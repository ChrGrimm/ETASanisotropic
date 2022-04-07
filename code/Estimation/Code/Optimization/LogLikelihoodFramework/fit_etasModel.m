function [ Catalog, ...
           IterationLog, ...
           paramETAS, ...
           SpatData, ...
           triggProbMatrix ] = fit_etasModel( Catalog, ...
                                               Inputs, ...
                                               ModelFuncs, ...
                                               SpatData, ...
                                               TempData, ...
                                               writeLog )
    
    %% Parameter Estimation for ETAS & ETAS-Incomplete models
    if ismember( Inputs.ModelSettings.modelType, {'ETAS', 'ETAS-Incomplete'} )
        %% Preparatory steps
        % Take square root of paramETAS (for numerical reasons)
        sqrtParamETAS = sqrt(Inputs.ModelSettings.iniParamETAS);
        % Initialize diagnostic variables
        [ IterationLog, hessMatrix, nIter ] = initialize_etasVariables( length(sqrtParamETAS) );    

        %% Loop over model iterations
        for iIter = 1:nIter
            
            if iIter==6
                test=1;
            end

            disp(['Iteration ', num2str(iIter), '/', num2str(nIter), ' ...'])        
            tStartIteration = tic;

            %% Declustering
            % Probabilistic declustering in ETAS
            [ Catalog, ...
              aggregBackgrRate ] = estimate_backgrProbabilities( Catalog, ...
                                                                  Inputs, ...
                                                                  ModelFuncs, ...
                                                                  sqrtParamETAS, ...
                                                                  SpatData, ...
                                                                  nIter-iIter+1, ...
                                                                  writeLog );

            %% Loglikelihood estimation of ETAS parameters
            [ IterationLog, ...
                hessMatrix, ...
                Catalog, ...
                sqrtParamETAS, ...
                SpatData, ...
                hasConverged ] = optimize_davidsonFletcherPowell( Catalog, ...
                                                                    Inputs, ...
                                                                    ModelFuncs, ...
                                                                    IterationLog, ...
                                                                    hessMatrix, ...
                                                                    sqrtParamETAS, ...
                                                                    aggregBackgrRate, ...
                                                                    SpatData, ...
                                                                    TempData, ...
                                                                    iIter, ...
                                                                    writeLog );

            tEndIteration = toc(tStartIteration)/60;
            disp(['Iteration ', num2str(iIter), ' took ', num2str(tEndIteration), ' minutes.'])

            % Stop estimation if current iteration ran into unrealistic parameter values
            if hasConverged
                break;           
            end
            
        end

        % Compute trigger probability matrix
        paramETAS               = sqrtParamETAS.^2;
        idxTarget               = find( Catalog.flag > 0 );
        ModelFuncs              = set_modelFunctions( Inputs, struct, sqrtParamETAS );
        [ ~,~, triggProbMatrix ]  = estimate_eventRates( Catalog, ...
                                                        Inputs, ...
                                                        sqrtParamETAS, ...
                                                        ModelFuncs, ...
                                                        idxTarget, ...
                                                        -1, ...
                                                        -1, -1, -1, -1, ...
                                                        false, ...
                                                        'time-space rates' );
        
        if strcmp( Inputs.ModelSettings.modelType, 'ETAS-Incomplete' )                                          
            F                               = estimate_spatialIntegral( 0, Catalog, Inputs, sqrtParamETAS(6:8).^2, ModelFuncs, SpatData );
            T2minusT1                       = Inputs.TargetSettings.tWindow_days(2)-Inputs.TargetSettings.tWindow_days(1); 
            nBackgrEventsPerDay             = aggregBackgrRate / T2minusT1;
            isTargetNonDupl                 = Catalog.flag > 0 & ~Catalog.isDupl;
            Catalog.N0_t(isTargetNonDupl)   = estimate_eventRates( Catalog, Inputs, sqrtParamETAS, ModelFuncs, idxTarget, ...
                                                                    nBackgrEventsPerDay, F, -1, -1, -1, false, 'time rates' );
        end
        
    end
    
    %% Magnitude Distribution Estimation for ETAS and Poisson models
    if ismember( Inputs.ModelSettings.modelType, {'ETAS', 'Poisson'} )
        % Fit of Gutenberg-Richter law
        paramETAS(10) = estimate_beta( Catalog, Inputs.TargetSettings.Mc );
    end

    %% Delete columns not needed to store
    Catalog = removevars(Catalog, {'isInExtent', 'dist', 'dist2'});
    
end
    