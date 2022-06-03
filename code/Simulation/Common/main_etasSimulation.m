function main_etasSimulation

    %% Seed for reproducability
    rng(1)
    
    %% Directory
    disp(['... Run ETAS Simulation of ', pwd, ' ...'])
       
    %% Load input data
    SimulInputs = load('Inputs_simulation');

    %% Call Simulation Functions
    if strcmp( SimulInputs.ExperimentSettings.TypeSimulation, 'Catalog' )
        % Synthetic Catalogs
        tAnalysis_synthCat = simulate_synthCatalogs( SimulInputs );
                                                        
    elseif strcmp( SimulInputs.ExperimentSettings.TypeSimulation, 'WPET' )
        % Synthetic Sequences
        simulate_synthWPET( SimulInputs );

    elseif strcmp( SimulInputs.ExperimentSettings.TypeSimulation, 'Sequence' )
        % Synthetic Sequences
        [ tAnalysis_synthSeq, ...
          tDetails_synthSeq ] = simulate_synthSequences( SimulInputs );
      
    else
        error('Unknown type of simulation')
        
    end
    
    %% Store results in .mat file
%     save Simulation_results.mat tAnalysis_synthCat 
%     cd('../Results')
%     filename = ['Simulation_results_', strNameSimulSettings, '.mat'];
%     save(filename, 'tAnalysis_synthCat')
                                                    
%     %% Plot results
%     plot_synthCatalog( SynthCatalogs, SimulationStatistic )    
