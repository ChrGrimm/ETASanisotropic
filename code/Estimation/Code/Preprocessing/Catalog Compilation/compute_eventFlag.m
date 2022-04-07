function eventFlag = compute_eventFlag( Catalog, TargetSettings )

    % boolean, if events are inside of polygon or on polygon edge
    inPoly = inpolygon( Catalog.x, Catalog.y, TargetSettings.polygonXY(:,1), TargetSettings.polygonXY(:,2) );
    
    %% Determine true/false, whether events occurred within target time window
    inTime = any( Catalog.t >= TargetSettings.tWindow_days(:,1)' & Catalog.t <= TargetSettings.tWindow_days(:,2)', 2 );
      
    %% Event flags
    eventFlag                        = zeros(size(Catalog,1),1);
    eventFlag( inPoly & inTime )     = 1;     % target event inside of polygon  
    eventFlag( ~inPoly & inTime )    = 0;     % event in target time but outside of polygon
    eventFlag( ~inPoly & ~inTime )   = -0.1;  % event before target time and outside of polygon
    eventFlag( inPoly & ~inTime )    = -1;    % event before target time but inside of polygon             
 
end