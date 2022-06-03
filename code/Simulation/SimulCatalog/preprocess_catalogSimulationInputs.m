function [ BackgrDistr, ...
           InitialEventSet, ...
           paramETAS, ...
           tAnalysis ] = preprocess_simulationInputs( SimulInputs, ...
                                                             InitialEventSet, ...
                                                             paramETAS, ...
                                                             ModelFuncs )

    %% Precompute background distribution
    BackgrDistr            = estimate_bkgdRatesOnGrid( TargetSettings, ExperimentSettings ); 
%     nBackgrEvents   = poissrnd( expectedValue, SimulSettings.nRealizations, 1 );
    
    %% Overwrite b value (magnitude distribution) 
    if isempty(ExperimentSettings.bValue_overwrite)
        disp(['Use estimated Gutenberg-Richter parameter \beta = b*log(10) = ', num2str(paramETAS(10))])
    elseif isnumeric(ExperimentSettings.bValue_overwrite)
        if length(ExperimentSettings.bValue_overwrite)==1
            paramETAS(10) = ExperimentSettings.bValue_overwrite*log(10);
            disp(['Manually set Gutenberg-Richter parameter to \beta = b*log(10) = ', num2str(paramETAS(10))])
        else
            ExperimentSettings.magnSample = ExperimentSettings.bValue_overwrite;
        end
    else
        error('Invalid choice of magnitude sampling method')
    end
    
    %% Compute branching ratio
    branchRatio = compute_branchingRatio( -1, -1, paramETAS(10), TargetSettings, ModelFuncs, 'simulation' );
    disp(['Simulation is performed with branching ratio = ', num2str(branchRatio)])
    
    %% Reformat InitialEventSet
    InitialEventSet(:,{'isDupl','date','lon','lat','tempRestr'}) = [];
    InitialEventSet.triggerID   = NaN(size(InitialEventSet,1),1);
    InitialEventSet.clusterID   = NaN(size(InitialEventSet,1),1);
    periodID                    = -1*ones(size(InitialEventSet,1),1);
    InitialEventSet             = addvars(InitialEventSet, periodID, 'Before', 'id');
    
    %% Initialize tAnalysis
    tAnalysis               = create_analysisTable4catalogs( ExperimentSettings.nRealizations );
                                                            
end