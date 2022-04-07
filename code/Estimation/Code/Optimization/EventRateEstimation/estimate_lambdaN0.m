function [ N0, ...
           dN0 ] = estimate_lambdaN0( Catalog, ...
                                       Inputs, ...
                                       nBackgrPerDay, ...
                                       sqrtParamETAS, ...
                                       idxTargetEvents, ...
                                       ModelFuncs, ...
                                       F, dF_D, dF_gamma, dF_q, ...
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
%   times2evaluate     : row indices of "Catalog" for which lambdas are supposed to be calculated

    %% Data Inputs and Initializations
    % Evaluate whether gradients should be computed
    computeGradients        = strcmp(strMode, 'with gradient');
    computeTriggerMatrix    = strcmp(strMode, 'trigger matrix');
    
    % Extract ETAS parameters (they are in square root format  for numerical reasons)
    [mu, A, ~, ~, p, ~, ~, ~, Tb] = extract_paramETASvalues( sqrtParamETAS, 'sqrt' );
    
    % Extract event-specific information in columns of table Catalog
    t           = Catalog.t;
    mag         = Catalog.mag;
    evWeight    = Catalog.evWeight;
    tempRestr   = Catalog.tempRestr;

    % Initialize output vectors and matrices
    [ ~, ~, ~, ...
        N0_trig, dN0 ] = initialize_eventRateOutputs( length(idxTargetEvents), ...
                                                      length(sqrtParamETAS), ...
                                                      size(Catalog,1), ...
                                                      computeGradients, ...
                                                      computeTriggerMatrix );
    
    %% Loop over events in times2evaluate
    Mc              = Inputs.TargetSettings.Mc;
    productivity    = ModelFuncs.fProductivity( mag, Mc );
    factor          = productivity .* evWeight .* F;
    factor_dD       = productivity .* evWeight .* dF_D;
    factor_dGamma   = productivity .* evWeight .* dF_gamma;
    factor_dQ       = productivity .* evWeight .* dF_q;
    
    fTempK_inn      = ModelFuncs.fTempK_inn;
    fixedParamETAS  = Inputs.ModelSettings.fixedParamETAS;
    spaceModel      = Inputs.ModelSettings.spaceModel;
    Mc              = Inputs.TargetSettings.Mc;
    
    for jCounter = 1:length(idxTargetEvents)
        
        isContribEv = t < t(idxTargetEvents(jCounter)) & ...
                      t >= t(idxTargetEvents(jCounter)) - tempRestr; % & t(idxTargetEvents(jCounter))-t <= 365;
        
        fTempK_inn_ij            = fTempK_inn(t(idxTargetEvents(jCounter))-t(isContribEv));
        fTempK_inn_ij_toMinusP   = fTempK_inn_ij.^(-p); 
                
        % Calculate trigger rates contributed from past events to event jCounter
        N0_trig_ij          = evWeight(jCounter) * factor(isContribEv) .* fTempK_inn_ij_toMinusP;
        N0_trig(jCounter)   = sum( N0_trig_ij );
        
        % Plausibility check
        if any(~isreal(N0_trig_ij))
            error('Imaginary result!')
        end
                           
        %% Computation of temporal and spatial event rate
        if computeGradients          
            
            if ~fixedParamETAS(2)
                dN0(2,jCounter)   = Tb * N0_trig(jCounter)/A;
            end
            if ~fixedParamETAS(3)
                dN0(3,jCounter)   = Tb * sum( N0_trig_ij .* (mag(isContribEv)-Mc) );
            end
            if ~fixedParamETAS(4)
                tmp_c           = -p ./ fTempK_inn_ij; % * ci(isContribEv)
                dN0(4,jCounter) = Tb * sum( N0_trig_ij .* tmp_c );
            end
            if ~fixedParamETAS(5)
                tmp_p           = -log(fTempK_inn_ij);
                dN0(5,jCounter) = Tb * sum( N0_trig_ij .* tmp_p );
            end
            if ~strcmp(spaceModel, 'none')
                if ~fixedParamETAS(6)
                    dN0(6,jCounter)  = Tb * sum( factor_dD(isContribEv) .* fTempK_inn_ij_toMinusP );
                end
                if ~fixedParamETAS(7)
                    dN0(7,jCounter)  = Tb * sum( factor_dGamma(isContribEv) .* fTempK_inn_ij_toMinusP );
                end
                if ~fixedParamETAS(8)
                    dN0(8,jCounter)  = Tb * sum( factor_dQ(isContribEv) .* fTempK_inn_ij_toMinusP );
                end
            end

        end
        
    end

    % Add time-homogeneous background rate
    N0          = Tb * ( mu * nBackgrPerDay + N0_trig );
    
    % Gradients
    if computeGradients
        if ~fixedParamETAS(1)
            if strcmp( spaceModel, 'none' )
                dN0(1,:)            = Tb;
            else
                dN0(1,:)            = Tb * nBackgrPerDay;
            end
        end
        if ~fixedParamETAS(9)
            dN0(9,:)    = N0/Tb;
        end
        dN0         = dN0 * 2 .* sqrtParamETAS;
    end
    
    %% Sum up N0 and derivatives for duplicated events   
    if any(Catalog.isDupl)
        
        idxDuplicates   = find( Catalog.id == Inputs.SpaceSettings.evIDs_twoStrikes );
        
        if ismember(strMode, {'w/o gradient', 'with gradient'})
            N0          = aggregate_duplicateEventRates( N0, idxDuplicates, 'col' );
            
            if computeGradients
                dN0    = aggregate_duplicateEventRates( dN0, idxDuplicates, 'col' );                
            end
            
        elseif computeTriggerMatrix
            N0_matrix       = aggregate_duplicateEventRates( N0_matrix, idxDuplicates, 'col' ); 
            N0_matrix       = aggregate_duplicateEventRates( N0_matrix, idxDuplicates, 'row' ); 
            
        end  
    end
    
end
