function Inputs = modify_inputs4clusterReEstimation( Catalog, Cluster, Inputs, paramETAS )

    %% Reduce catalog to cluster events
    Catalog = Catalog( ismember(Catalog.id, Cluster.eventIDs{1}) & ~Catalog.isDupl, : );
    save Cluster_Catalog.mat Catalog
%     save( ['Cat_Cluster', num2str(iCluster), '.mat'], 'Catalog' )
     
    %% Create circular polygon around mainshock, including all cluster events
    isAniso         = Catalog.mag >= Inputs.SpaceSettings.anisoFromMw;
    maxDistance     = max( getDistance( Cluster.mainshLon, Cluster.mainshLat, Catalog.lon, Catalog.lat, 'km' ) + ...
                           isAniso .* estimate_rupSize( Catalog.mag, Inputs.SpaceSettings.tectonicType, Inputs.SpaceSettings.faultingStyle, 'km' ) );
    radiusPolygon   = ceil(maxDistance) + 1;
    polygonTarget   = create_circularPolygon( radiusPolygon, Cluster.mainshLon, Cluster.mainshLat );
    polygonComple   = polygonTarget;
    save Cluster_Polygon.mat polygonTarget polygonComple
    
    %% Modify inputs
    % Model settings
    ModelSettings               = Inputs.ModelSettings;
    ModelSettings.mu_ini        = paramETAS(1);
    ModelSettings.A_ini         = paramETAS(2);
    ModelSettings.alpha_ini     = paramETAS(3);
    ModelSettings.c_ini         = paramETAS(4);
    ModelSettings.p_ini         = paramETAS(5);
    ModelSettings.D_ini         = paramETAS(6)*111^2;
    ModelSettings.gamma_ini     = paramETAS(7);
    ModelSettings.q_ini         = paramETAS(8);
    ModelSettings.Tb_ini        = paramETAS(9)*24*60*60;
    ModelSettings.b_ini         = paramETAS(10)/log(10);
    ModelSettings.iniBackgrProb = Catalog.backgrProb;
    
    % Space settings (no changes)
    SpaceSettings               = Inputs.SpaceSettings;
    
    % Target settings (modify catalog, polygon and time window)
    TargetSettings              = Inputs.TargetSettings;
    TargetSettings.strCatalog   = 'Cluster_Catalog.mat';
    TargetSettings.strPolygon   = 'Cluster_Polygon.mat';
%     TargetSettings.tIni         = TargetSettings.tStart; 
    TargetSettings.tIni         = Catalog.date(1) - days(10);
    TargetSettings.tStart       = Catalog.date(1) - days(10);
    TargetSettings.tEnd         = Catalog.date(end);
    
    % Time settings (no changes)
    TimeSettings                = Inputs.TimeSettings;
    
    save Inputs_estimation.mat ModelSettings SpaceSettings TargetSettings TimeSettings

end