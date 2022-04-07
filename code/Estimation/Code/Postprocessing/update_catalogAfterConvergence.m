function Catalog = update_catalogAfterConvergence( Catalog, Inputs, ModelFuncs, paramETAS )
    
    %% Update event rate and background probability after convergence
    isTarget                    = Catalog.flag  > 0;
    Catalog.lambda(isTarget)    = estimate_lambda( Catalog, ...
                                                   Inputs, ...
                                                   sqrt(paramETAS), ...
                                                   find(isTarget), ...
                                                   ModelFuncs, ...
                                                   'decluster' );
                                               
    Catalog.backgrProb(isTarget) = paramETAS(1) * Catalog.backgrRate(isTarget) ./ Catalog.lambda(isTarget);
    
end

%% Old code pieces    
%     %% Compute integral of spatial trigger function per event
%     % Extract spatial parameters
%     spatParamETAS = sqrtParamETAS(6:8).^2;
%     
%     % Estimate integral per event
%     Catalog.spatIntegr   = estimate_spatialIntegral(0, ...
%                                                     Catalog, ...
%                                                     Inputs, ...
%                                                     spatParamETAS, ...
%                                                     ModelFuncs, ...
%                                                     SpatData);
    
%     %% Compute bandwidth and background integral per event    
%     % Initialize new columns
%     Catalog.backgrIntegr            = nan(size(Catalog.id));    
%     isTargetT                       = Catalog.flag >= 0;
%     Catalog.backgrIntegr(isTargetT) = SpatData.backgrIntegral;
                                                            
%     % Set event and background rates of non-target events to NaN
%     Catalog.lambda( ~isTarget )                 = NaN;
%     Catalog.backgrRate( ~isTarget )             = NaN;
