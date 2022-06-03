function [ Inputs, ...
           BackgrDistr, ...
           EventListWPET, ...
           InputsWPET ] = preprocess_simulationInputs( Inputs )
        
    %% Overwrite b value (magnitude distribution)
    if isscalar(Inputs.ExperimentSettings.sampleMagnitudes)
        Inputs.paramETAS(10) = log(10) * Inputs.ExperimentSettings.sampleMagnitudes;
        disp(['Manually set Gutenberg-Richter parameter to \beta = b*log(10) = ', num2str(Inputs.paramETAS(10))])
    end
    
    %% Precompute background distribution (catalog simulation only)
    if strcmp( Inputs.ExperimentSettings.TypeSimulation, 'Catalog' )
        BackgrDistr = estimate_bkgdRatesOnGrid( Inputs.TargetSettings, ...
                                                Inputs.ExperimentSettings ); 
    else
        BackgrDistr = 'dummy';
    end
    
    %% Compute and display branching ratio
    if isempty(Inputs.ExperimentSettings.sampleMagnitudes) || isscalar(Inputs.ExperimentSettings.sampleMagnitudes) 
        branchRatio = compute_branchingRatio( -1, -1, ...
                                              Inputs.paramETAS(10), ...
                                              Inputs.TargetSettings, ...
                                              Inputs.ModelFuncs, ...
                                              'simulation' );
        disp(['Simulation is performed with branching ratio = ', num2str(branchRatio)])
    end
   
    %% Compile initial event set
    if strcmp( Inputs.ExperimentSettings.TypeSimulation, 'WPET' )
        InputsWPET = load('Inputs_wpet.mat', 'ELT', 'MET', 'Faults', 'Zones');
        load('Inputs_wpet.mat', 'EventListWPET');
%         if strcmp(Inputs.ExperimentSettings.sampleMagnitudes, 'MET')
%             InputsWPET.isActiveFault = evaluate_activeFaultsPerPeriod( InputsWPET.MET, ...
%                                                                        InputsWPET.Faults, ...
%                                                                        EventListWPET );
%         end
    else
        EventListWPET   = 'dummy';
        InputsWPET      = 'dummy';
    end
    
end