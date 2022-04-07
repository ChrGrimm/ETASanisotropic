function Catalog = delete_eventsOutsidePolygon( Catalog, Inputs )

    %% Initializations
    % Extract settings
    polygon         = Inputs.TargetSettings.polygonXY;
    spaceUnit       = Inputs.SpaceSettings.spaceUnit;
    
    % Initialize variables
    idxEvents       = find( abs(Catalog.flag)<0.5 );
    eventCoords     = [Catalog.x(idxEvents)'; Catalog.y(idxEvents)'];
    nPolyEdges      = size(polygon,1)-1;
    dists2polygon   = zeros( length(idxEvents), nPolyEdges );

    %% Compute distances between events and polygon edges
    for iEdge = 1:nPolyEdges
        
        dists2polygon(:,iEdge) = get_dist2line(polygon(iEdge,:), polygon(iEdge+1,:), eventCoords, spaceUnit, 0);
        
    end
    
    %% Delete events, whose closest distance to polygon exceeds spatial extent
    idx2delete                  = idxEvents( min(dists2polygon, [], 2) > Catalog.spatRestr(idxEvents) );
    Catalog( idx2delete, : )    = [];
    %(  & strcmp(Catalog.typeKernel, 'full'), : ) = [];
    