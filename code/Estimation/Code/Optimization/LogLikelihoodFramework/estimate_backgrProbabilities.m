function [ Catalog, ...
            aggregBackgrRate ] = estimate_backgrProbabilities( Catalog, ...
                                                              Inputs, ...
                                                              ModelFuncs, ...
                                                              sqrtParamETAS, ...
                                                              SpatData, ...
                                                              nIterations, ...
                                                              writeLog )
      
    if strcmp(Inputs.ModelSettings.spaceModel, 'none')
        %% If no space model
        Catalog.backgrRate(:)   = Catalog.flag == 1;
        aggregBackgrRate        = sum( Inputs.TargetSettings.tWindow_days(:,2)-Inputs.TargetSettings.tWindow_days(:,1) );
            
    else
        %% Default: With space model
        timeWindow_days = Inputs.TargetSettings.tWindow_days;
        
        for iDecl = 1:nIterations
            if writeLog
                % Print iteration number
                disp(['Declustering ', num2str(iDecl), '/', num2str(nIterations), ' ...'])                                                   
            end
            
            % Background rate per event
            isTargetT                       = Catalog.flag >= 0 & ~Catalog.isDupl;
            isTarget                        = Catalog.flag > 0 & ~Catalog.isDupl;
            Catalog.backgrRate(isTarget)    = sum( Catalog.backgrProb(isTargetT) .* SpatData.gaussDensity )' ...
                                                    / sum(timeWindow_days(:,2)-timeWindow_days(:,1)); % / diff(timeWindow_days);
            if any(Catalog.isDupl)
                idxDuplicate                    = find( Catalog.isDupl );
                Catalog.backgrRate(idxDuplicate)= Catalog.backgrRate( idxDuplicate-1 );
                Catalog.backgrRate              = Catalog.backgrRate .* Catalog.evWeight;
            end

            % Last iteration: Compute integrated background rate over target time and space
            if iDecl == nIterations
                aggregBackgrRate = sum( sum( Catalog.backgrProb(isTargetT) .* SpatData.backgrIntegral ) );
            else
                aggregBackgrRate = -1; % dummy value
            end

            %% 1st iteration: Compute pure trigger rates lambda
            isTarget = Catalog.flag > 0;
            if iDecl == 1
                % Define functions with current parameters
                ModelFuncs = set_modelFunctions( Inputs, ModelFuncs, sqrtParamETAS );
                
                % Compute compute lambdas                                           
                lambda = estimate_eventRates( Catalog, ...
                                              Inputs, ...
                                              sqrtParamETAS, ...
                                              ModelFuncs, ...
                                              find(isTarget), ...
                                              -1, ...
                                              -1, -1, -1, -1, ...
                                              false, ...
                                              'decluster' );
                                       
                % Subtract background rate to obtain pure trigger rate
                sum_lambda_trig = lambda' - sqrtParamETAS(1)^2 * Catalog.backgrRate(isTarget); 
            end

            lambda                          = sum_lambda_trig + sqrtParamETAS(1)^2 * Catalog.backgrRate(isTarget);
            Catalog.backgrProb(isTarget)    = sqrtParamETAS(1)^2 * Catalog.backgrRate(isTarget) ./ lambda;
            Catalog.lambda(isTarget)        = lambda;
            
        end
    end

end

% Old code pieces
%     Catalog.backgrRate(isTarget)    = sum( Catalog.backgrProb(isTargetT) .* SpatData.gaussDensity )' ...
%                                     / sum(timeWindow_days(:,2)-timeWindow_days(:,1)) / sqrtParamETAS(1)^2;
%   Catalog.backgrRate    = sum( Catalog.backgrProb .* SpatData.gaussDensity )' / diff(timeWindow);

