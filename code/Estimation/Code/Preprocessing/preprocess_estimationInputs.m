function [ Catalog, ...
            Inputs, ...
            ModelFuncs, ...
            SpatData, ...
            TempData ] = preprocess_estimationInputs( Inputs, writeLog )
        
    %% Load original catalog
    load( Inputs.TargetSettings.pathCatalog, 'Catalog' );
    
    %% Modify and check ETAS model inputs
    [Inputs, ModelFuncs] = preprocess_estimationSettings( Inputs );
                                                       
    % Print inputs
    if writeLog
        print_estimationSettings( Inputs );
    end
    
    %% Compile ETAS input catalog
    Catalog = compile_catalog4etas( Inputs, ModelFuncs, Catalog, 'estimation' ); 
    
    % Precompute spatial background distribution & polygon grids
    [ SpatData, Catalog ] = precompute_spatialData( Catalog, Inputs, ModelFuncs );
                                                
    % Precompute evaluation times for N0 estimation (ETAS-incomplete only!)
    TempData = precompute_temporalData( Catalog, ...
                                        Inputs, ...
                                        ModelFuncs, ...
                                        SpatData );
    
end