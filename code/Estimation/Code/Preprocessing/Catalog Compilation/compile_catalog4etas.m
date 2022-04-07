function Catalog = compile_catalog4etas( Inputs, ModelFuncs, origCatalog, version )

    %% Filter original catalog
    origCatalog     = filter_origCatalog( origCatalog, ...
                                          Inputs.TargetSettings, ...
                                          version );
    nEvents         = size(origCatalog, 1);
                                                         
    %% Initialize ETAS Catalog
    Catalog         = table;
    Catalog.id      = origCatalog.id;
    Catalog.isDupl  = false(nEvents,1);
    
    % Time Info
    Catalog.date    = origCatalog.date;
    Catalog.t       = days(origCatalog.date - Inputs.TargetSettings.tWindow(1,1));
    
    % Space Info
    Catalog.lon     = origCatalog.lon;
    Catalog.lat     = origCatalog.lat;
    Catalog.x       = ModelFuncs.fConvert_lon2x( Catalog.lon );
    Catalog.y       = ModelFuncs.fConvert_lat2y( Catalog.lat );
    
    % Depth and Magnitude
%     Catalog.depth   = origCatalog.depth;
    Catalog.mag     = origCatalog.mag;  
    
    % Event flag and weight
    Catalog.flag            = compute_eventFlag( Catalog, Inputs.TargetSettings );
    Catalog.evWeight        = ones(nEvents,1);
    if strcmp(version, 'estimation')
        Catalog.lambda          = zeros(nEvents, 1);
        Catalog.backgrRate      = zeros(nEvents, 1);
        Catalog.backgrProb(:)   = Inputs.ModelSettings.iniBackgrProb;
    end
    
    % Spatial kernel characteristics for isotropic models
    Catalog.typeKernel(:)   = {'iso'};
    Catalog.rupExtent       = zeros(nEvents,1);
    Catalog.tempRestr       = ModelFuncs.fRestrTime_days( Catalog.mag );
    Catalog.spatRestr       = ModelFuncs.fRestrSpace_degrees( Catalog.mag, 'iso' );
    Catalog.strike          = nan(nEvents,1);
    Catalog.epiPos          = nan(nEvents,1);
    
    % Spatial kernel characteristics for anisotropic models
    if strcmp(Inputs.ModelSettings.spaceModel, 'aniso')
        % Overwrite information for anisotropic trigger events
        idxAniso                        = find( Catalog.mag >= Inputs.SpaceSettings.anisoFromMw );
        Catalog.typeKernel(idxAniso)    = {'aniso'};
        Catalog.rupExtent(idxAniso)     = ModelFuncs.fRupLength( Catalog.mag(idxAniso) );
        Catalog.spatRestr(idxAniso)     = ModelFuncs.fRestrSpace_degrees( Catalog.mag(idxAniso), 'aniso' );
                                                                     
        % Create event duplicate, if two strike directions should be modeled
        Catalog = manage_eventsWithTwoStrikes( Catalog, Inputs.SpaceSettings.evIDs_twoStrikes, version );
        
        % Update indices of anisotropic events
        idxAniso = find( Catalog.mag >= Inputs.SpaceSettings.anisoFromMw );
        
        % Rupture orientation
        if ~isempty(idxAniso)
            Catalog = estimate_nodalPlaneSolutions( Catalog, origCatalog, ...
                                                    Inputs, ModelFuncs, idxAniso );
        end
    end
   
    %% Delete outside events whose spatial extent does not reach modeled spatial window
    if strcmp(version, 'estimation')
        Catalog = delete_eventsOutsidePolygon( Catalog, Inputs );
    end
    
    %% Cross-event distances and relations 
    if strcmp(version, 'estimation')
        % - boolean that marks which events are within spatio-temporal extent
        % - Number of events that are within spatio-temporal extent
        % - (Squared) distances to all events within spatio-temporal extent
        % - Start and end point of rupture line (only needed for aniso model; dummy otherwise)
        [ Catalog.isInExtent, ...
          Catalog.dist, ...
          Catalog.dist2, ...      
          Catalog.rupStart, ...
          Catalog.rupEnd ] = precompute_crossEventDistances( Catalog, ...
                                                             Inputs.SpaceSettings.spaceUnit, ...
                                                             Inputs.SpaceSettings.minXeventDist );
    end
    
    %% Check that no anisotropic rupture line crosses polygon boundaries
    for idxAniso = find(strcmp(Catalog.typeKernel', 'aniso'))
        intersecX = polyxpoly( Inputs.TargetSettings.polygonXY(:,1), ...
                                Inputs.TargetSettings.polygonXY(:,2), ...
                                [Catalog.rupStart(idxAniso,1), Catalog.rupEnd(idxAniso,1)], ...
                                [Catalog.rupStart(idxAniso,2), Catalog.rupEnd(idxAniso,2)] );
        if ~isempty(intersecX)    
            % Send warning message
            disp(['The rupture segment of event ', num2str(Catalog.id(idxAniso)), ' was partly inside and outside of target polygon. ', ...
                     'Polygon was modified to cover entire rupture segment'])
                 
            % Create very narrow polygon around rupture segment
            pointsXY                        = [Catalog.rupStart(idxAniso,:); Catalog.rupEnd(idxAniso,:)];
            polyAroundRupSegm               = polybuffer(pointsXY,'lines',10^-5);
            
            unionPolygon                    = union(polyshape(Inputs.TargetSettings.polygonXY), polyAroundRupSegm);
            Inputs.TargetSettings.polygonXY = unionPolygon.Vertices;
            
            if Catalog.flag(idxAniso) == 0
                Catalog.flag(idxAniso) = 1;
            elseif Catalog.flag(idxAniso) == -0.1
                Catalog.flag(idxAniso) = -1;
            end
            
%             % Plot situation
%             figure
%             plot(Inputs.TargetSettings.polygonXY(:,1), Inputs.TargetSettings.polygonXY(:,2))
%             hold on
%             plot( [Catalog.rupStart(idxAniso,1), Catalog.rupEnd(idxAniso,1)], ...
%                   [Catalog.rupStart(idxAniso,2), Catalog.rupEnd(idxAniso,2)] )
            
        end
    end
  
end
