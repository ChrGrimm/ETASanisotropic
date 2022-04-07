function [ R0, ...
           triggProbMatrix, ...
           dR0, ...
           N0, ...
           dN0 ] = estimate_lambda( Catalog, ...
                                    Inputs, ...
                                    sqrtParamETAS, ...
                                    idxTargetEvents, ...
                                    ModelFuncs, ...
                                    strMode )
% 
% This function calculates the ETAS rate lambda at the exact historical
% event dates and locations, either with isotropic or anisotropic spatial
% distribution, optionally including calculation of gradiants.
%
% INPUT: 
%   inclGradients    : boolean, states whether gradiants are to be calculated as well
%   Catalog     : table of historical event data
%                   .t: event time in days after initial event
%                   .x: centralized, scaled longitude coordinate
%                   .y: centralized latitude coordinate
%   modelDesign:    string      
%   sqrtParamETAS   : vector of current parameter estimates: (mu, A, c, alpha, p, D, q, gamma)
%   idxTargetEvents     : row indices of "Catalog" for which lambdas are supposed to be calculated

    %% Data Inputs and Initializations
    % Evaluate whether gradients should be computed
    computeGradients        = strcmp(strMode, 'with gradient');
    computeTriggerMatrix    = strcmp(strMode, 'trigger matrix');
    
    % Extract ETAS parameters (they are in square root format  for numerical reasons)
    [mu, A, ~, ~, p, D, ~, q] = extract_paramETASvalues( sqrtParamETAS, 'sqrt' );
    
    % Extract event-specific information in columns of table Catalog
    t           = Catalog.t;
    mag         = Catalog.mag;
    isInExtent  = Catalog.isInExtent;
    dist        = Catalog.dist;
    dist2       = Catalog.dist2;
    rupExtent   = Catalog.rupExtent;
    spatRestr   = Catalog.spatRestr;
    tempRestr   = Catalog.tempRestr;
    evWeight    = Catalog.evWeight;
    
    Mc          = Inputs.TargetSettings.Mc;
    spaceModel  = Inputs.ModelSettings.spaceModel;

    % Initialize output vectors and matrices
    [ R0, dR0, R0_matrix, N0, dN0 ] = initialize_eventRateOutputs( length(idxTargetEvents), ...
                                                                  length(sqrtParamETAS), ...
                                                                  size(Catalog,1), ...
                                                                  computeGradients, ...
                                                                  computeTriggerMatrix );
    triggProbMatrix = -1;
    
    % Compute normalization of spatial kernel
    fSpatK_inn  = ModelFuncs.fSpatK_inn;
    fTempK_inn  = ModelFuncs.fTempK_inn;
    factor      = ModelFuncs.fProductivity( mag, Mc ) ...
                    .* ModelFuncs.fSpatK_factor( spatRestr, mag, rupExtent ) ...
                    .* ModelFuncs.fTempK_factor( tempRestr );
    fixedParamETAS = Inputs.ModelSettings.fixedParamETAS;
    
    [ norm_D, norm_q ] = compute_normaliz4evRateDeriv( fSpatK_inn, spatRestr, mag, rupExtent, q );    

    if computeGradients
        if ~fixedParamETAS(1)
            dR0(1,:) = Catalog.backgrRate(idxTargetEvents)';
        end
%         spatDistrWidth      = D*exp(gamma*(Catalog.mag-Mc));
%         isWiderThanMinWidth = spatDistrWidth > SpatKernel.minKernelWidth;
    end
    
    %% Loop over events in idxTargetEvents
    jCounter = 0;
    for jTargetEv = idxTargetEvents'
        % Update event counter
        jCounter = jCounter + 1;
        % Extract boolean marking events that contribute trigger rate
        isContribEv = isInExtent{jTargetEv};
        
        fTempK_inn_ij    = fTempK_inn(t(jTargetEv)-t(isContribEv));
        fSpatK_inn_ij    = fSpatK_inn(dist{jTargetEv}, dist2{jTargetEv}, mag(isContribEv), rupExtent(isContribEv));
        
        % Calculate trigger rates contributed from past events to event jCounter
        % Event weights of triggering and triggered event times incompleteness weights
        weights     = evWeight(isContribEv) * evWeight(jTargetEv);
        R0_trig     = weights .* factor(isContribEv) .* fTempK_inn_ij.^(-p) .* fSpatK_inn_ij.^(-q);
        
        % Plausibility check
        if any(~isreal(R0_trig))
            error('Imaginary result!')
        end
        
        % Sum up trigger rates from all past events, contributing to the event rate of event jCounter
        R0(jCounter) = sum( R0_trig );
        
        if computeTriggerMatrix && ~isempty(R0_trig)        
            % Store detailed trigger rates in matrix  
            R0_matrix(isContribEv,jCounter)    = R0_trig;
        end
                   
        %% Computation of temporal and spatial event rate
        if computeGradients
            
            % Gradiants by productivity parameters
            if ~fixedParamETAS(2)
                dR0(2,jCounter)   = R0(jCounter)/A;
            end
            if ~fixedParamETAS(3)
                dR0(3,jCounter)   = sum( R0_trig .* (mag(isContribEv)-Mc) );
            end
            if ~fixedParamETAS(4)
                tmp_c           = -p ./ fTempK_inn_ij;
                dR0(4,jCounter) = sum( R0_trig .* tmp_c );
            end
            if ~fixedParamETAS(5)
                tmp_p           = -log(fTempK_inn_ij);
                dR0(5,jCounter) = sum( R0_trig .* tmp_p );
            end
            if ~strcmp(spaceModel, 'none')
                % norm_D = norm_q = 0 if spatial kernel no normalized (see before start of loop)          
                h_D              = (q * (1-fSpatK_inn_ij.^(-1)) - 1 - norm_D(isContribEv)); % .* isWiderThanMinWidth(isContribEv);
                if ~fixedParamETAS(6)                
                    dR0(6,jCounter)  = sum( R0_trig .* 1/D .* h_D );
                end
                if ~fixedParamETAS(7)
                    dR0(7,jCounter)  = sum( R0_trig .* (mag(isContribEv)-Mc) .* h_D );
                end
                if ~fixedParamETAS(8)
                    h_q                         = (1/(q-1) - log(fSpatK_inn_ij) - norm_q(isContribEv));
                    dR0(8,jCounter)  = sum( R0_trig .* h_q );
                end
            end
        end
        
    end
    
    % Add time-homogeneous background rate
    R0     = R0 + mu * Catalog.backgrRate(idxTargetEvents)';
    dR0    = dR0 * 2 .* sqrtParamETAS;
    
    %% Sum up R0 and derivatives for duplicated events  
    if any(Catalog.isDupl)
        
        idxDuplicates   = find( Catalog.id == Inputs.SpaceSettings.evIDs_twoStrikes );
        
        if ismember(strMode, {'w/o gradient', 'with gradient'})
            R0              = aggregate_duplicateEventRates( R0, idxDuplicates, 'col' );
            if computeGradients
                dR0         = aggregate_duplicateEventRates( dR0, idxDuplicates, 'col' );                
            end
            
        elseif computeTriggerMatrix
            R0              = aggregate_duplicateEventRates( R0, idxDuplicates, 'col' );
            R0_matrix       = aggregate_duplicateEventRates( R0_matrix, idxDuplicates, 'col' ); 
            R0_matrix       = aggregate_duplicateEventRates( R0_matrix, idxDuplicates, 'row' );
            triggProbMatrix = R0_matrix ./ R0;
            
        end  
    end
    
end
