function print_estimationSettings( Inputs )
     
    %% Model Settings
    disp(' ETAS INPUTS')
    disp(' ... ModelSettings ... ')                            
    disp(['modelType        = ', Inputs.ModelSettings.modelType])                            
    disp(['spaceModel       = ', Inputs.ModelSettings.spaceModel])
    disp(['iniParamETAS     = ', num2str(Inputs.ModelSettings.iniParamETAS')]) 
    
    %% Target Settings
    disp(' ... TargetSettings ... ') 
    disp(['pathCatalog   = ', Inputs.TargetSettings.pathCatalog])
    disp(['pathPolygon   = ', Inputs.TargetSettings.pathPolygon])
    disp(['timeWindow from ', datestr(Inputs.TargetSettings.tWindow(1)), ' till ', ...
                              datestr(Inputs.TargetSettings.tWindow(2))])
    disp(['preHistory from ', datestr(Inputs.TargetSettings.tIni)])
    disp(['Mc            = ', num2str(Inputs.TargetSettings.Mc)])
    disp(['Mmax          = ', num2str(Inputs.TargetSettings.Mmax)])
    disp(['maxDepth      = ', num2str(Inputs.TargetSettings.maxDepth)])
    
    %% Spatial Kernel
    disp(' ... SpaceSettings ... ') 
    disp(['restrFactor           = ', num2str(Inputs.SpaceSettings.restrFactor)])
    disp(['minRestr              = ', num2str(Inputs.SpaceSettings.minRestr)])

    if strcmp( Inputs.ModelSettings.spaceModel, 'aniso' )
        disp(['anisoFromMw           = ', num2str(Inputs.SpaceSettings.anisoFromMw)])
        disp(['tStrikeEstim_hrs      = ', num2str(Inputs.SpaceSettings.tStrikeEstim_hrs)])
        disp(['evIDs_twoStrikes      = ', num2str(Inputs.SpaceSettings.evIDs_twoStrikes)])
        disp(['tectonicType          = ', Inputs.SpaceSettings.tectonicType])
        disp(['faultingStyle         = ', Inputs.SpaceSettings.faultingStyle])
        
    end
    
    %% Omori Law
    disp(' ... TimeSettings ... ') 
    disp(['restr_days     = ', num2str(Inputs.TimeSettings.restr_days)])
    
end