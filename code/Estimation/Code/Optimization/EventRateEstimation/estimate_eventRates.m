function [ eventRates, ...
           eventRates_grad, ...
           triggProbMatrix ] = estimate_eventRates( Catalog, ...
                                                    Inputs, ...
                                                    sqrtParamETAS, ...
                                                    ModelFuncs, ...
                                                    idxTargetEvents, ...
                                                    nBackgrPerDay, ...
                                                    F, dF_D, dF_gamma, dF_q, ...
                                                    computeGradients, ...
                                                    strModus )
    
    %% Initialize Outputs
    nTargetEvents   = length(idxTargetEvents);
    eventRates      = zeros( 1, nTargetEvents );
    eventRates_grad = zeros( length(sqrtParamETAS), nTargetEvents);
    triggRateMatrix = zeros( size(Catalog,1), nTargetEvents );
    
    %% Extract data to avoid time-consuming data processing in event loop
    % Extract ETAS parameters
    [mu, A, ~, ~, p, D, ~, q, Tb] = extract_paramETASvalues( sqrtParamETAS, 'sqrt' );
    
    % Extract settings
    Mc              = Inputs.TargetSettings.Mc;
    fixedParamETAS  = Inputs.ModelSettings.fixedParamETAS;
    spaceModel      = Inputs.ModelSettings.spaceModel;
    
    % Extract event-specific information in columns of table Catalog
    t           = Catalog.t;
    mag         = Catalog.mag;
    tempRestr   = Catalog.tempRestr;
    evWeight    = Catalog.evWeight;
    if ismember(strModus, {'time-space rates', 'decluster'})
        % in time-space version only
        isInExtent  = Catalog.isInExtent;
        dist        = Catalog.dist;
        dist2       = Catalog.dist2;
        rupExtent   = Catalog.rupExtent;
        backgrRate  = Catalog.backgrRate;
        spatRestr   = Catalog.spatRestr;
    end
    
    % Extract model functions
    fTempK_inn      = ModelFuncs.fTempK_inn;
    fSpatK_inn      = ModelFuncs.fSpatK_inn;
    productivity    = ModelFuncs.fProductivity( mag, Mc );
    
    %% Precomputations
    if ismember(strModus, {'time-space rates', 'decluster'})
        factor      = productivity ...
                        .* ModelFuncs.fSpatK_factor( spatRestr, mag, rupExtent ) ...
                        .* ModelFuncs.fTempK_factor( tempRestr );
        [ norm_D, ...
          norm_q ]  = compute_normaliz4evRateDeriv( fSpatK_inn, spatRestr, mag, rupExtent, q );    
                    
    elseif strcmp(strModus, 'time rates')
        factor          = productivity .* F;
        factor_dD       = productivity .* evWeight .* dF_D;
        factor_dGamma   = productivity .* evWeight .* dF_gamma;
        factor_dQ       = productivity .* evWeight .* dF_q;
        
    end
    
    %% EVENT RATE COMPUTATION (starts here)
    % Loop over events in idxTargetEvents
    for iEvent = 1:length(idxTargetEvents)
        
        iTargetEvent    = idxTargetEvents(iEvent);
        
        if ismember(strModus, {'time-space rates', 'decluster'})
            [ eventRates(iEvent), ...
              eventRates_grad(:,iEvent), ...
              triggRateMatrix(:,iEvent) ] = compute_singleEventRate( iTargetEvent, ...
                                                                     t, mag, tempRestr, evWeight, isInExtent, dist, dist2, rupExtent, backgrRate, ...
                                                                     fSpatK_inn, fTempK_inn, ...
                                                                     mu, A, p, D, q, Tb, ...
                                                                     Mc, -1, spaceModel, ...
                                                                     norm_D, norm_q, ...
                                                                     factor, -1, -1, -1, ...
                                                                     computeGradients, fixedParamETAS, ...
                                                                     strModus );
                                                                 
        else
            [ eventRates(iEvent), ...
              eventRates_grad(:,iEvent), ...
              triggRateMatrix(:,iEvent) ] = compute_singleEventRate( iTargetEvent, ...
                                                                     t, mag, tempRestr, evWeight, -1, -1, -1, -1, -1, ...
                                                                     fSpatK_inn, fTempK_inn, ...
                                                                     mu, A, p, D, q, Tb, ...
                                                                     Mc, nBackgrPerDay, spaceModel, ...
                                                                     -1, -1, ...
                                                                     factor, factor_dD, factor_dGamma, factor_dQ, ...
                                                                     computeGradients, fixedParamETAS, ...
                                                                     strModus );
        end
                                                                 
    end
    
            
    %% Inner derivative
%     if strcmp(strModus, 'time-space rates')
        eventRates_grad = eventRates_grad * 2 .* sqrtParamETAS;
%     end
    
    %% Sum up values for duplicated events  
    if any(Catalog.isDupl & Catalog.flag>0) && ~strcmp(strModus, 'decluster')
        idxDuplicates   = find( Catalog.id(Catalog.flag>0) == Inputs.SpaceSettings.evIDs_twoStrikes );
        eventRates      = aggregate_duplicateEventRates( eventRates, idxDuplicates, 'col' );
        eventRates_grad = aggregate_duplicateEventRates( eventRates_grad, idxDuplicates, 'col' );
        triggRateMatrix = aggregate_duplicateEventRates( triggRateMatrix, idxDuplicates, 'col' ); 
        triggRateMatrix = aggregate_duplicateEventRates( triggRateMatrix, idxDuplicates, 'row' );
    end
    
    %% Compute trigger probability matrix
    triggProbMatrix = triggRateMatrix ./ eventRates;                
                                                
end