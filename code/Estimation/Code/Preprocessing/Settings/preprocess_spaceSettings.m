function [ SpaceSettings, ...
            ModelFuncs ] = preprocess_spaceSettings( SpaceSettings, ...
                                                     spaceModel, ...
                                                     ModelFuncs )

    %% Modify space settings
    % Set space unit (km not implemented)
    SpaceSettings.spaceUnit         = 'degree';
    SpaceSettings.minRestr          = 0.2;
    SpaceSettings.minXeventDist     = 0.05;
    SpaceSettings.detailedIntegral  = true;
    SpaceSettings.sampleStrikes     = 0:1:179;
    SpaceSettings.sampleEpicenter   = 0:0.01:1;
    
    % Adjust inputs to space unit
    if strcmp(SpaceSettings.spaceUnit, 'degree')
        SpaceSettings.minRestr       = SpaceSettings.minRestr / 111;
        SpaceSettings.minXeventDist  = SpaceSettings.minXeventDist / 111;
    end
    
    %% Set model functions
    % Spatial restriction
    % Rupture length and anisotropic extension
    ModelFuncs.fRupLength          = @(mag) estimate_rupSize( mag, SpaceSettings.tectonicType, SpaceSettings.faultingStyle, 'degree' );
    % Spatial restriction
    ModelFuncs.fRestrSpace_factor  = @(typeKernel) SpaceSettings.restrFactor - 0.5*strcmp(typeKernel,'aniso');
    ModelFuncs.fRestrSpace_degrees = @(mag,typeKernel) max( SpaceSettings.minRestr, ModelFuncs.fRestrSpace_factor(typeKernel) .* ModelFuncs.fRupLength(mag) );

    %% Set minimum kernel width
    SpaceSettings.minKernelWidth = 0;
    
    % Fix settings for time-only models
    if strcmp( spaceModel, 'none' ) 
        SpaceSettings.detailedIntegral   = false;
        SpaceSettings.restrFactor        = 10^60;
    end
    
    if ~strcmp( spaceModel, 'aniso' )
        SpaceSettings.anisoFromMw = 99;
    end

end