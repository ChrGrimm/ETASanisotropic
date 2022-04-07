function Catalog = filter_origCatalog( origCatalog, TargetSettings, version )
             
    %% Preparation Steps
    % Round magnitudes and sort chronologically
    Catalog     = sortrows(origCatalog, 'date');
    Catalog.mag = round(Catalog.mag,1);
    
    % Remove duplicate events
    Catalog     = remove_eventDuplicates( Catalog );

    %% Filter by event times
    % Cut out events outside of time window
    if strcmp(version, 'estimation') || strcmp(version, 'benchmark')
        timeWindow  = [ TargetSettings.tIni, TargetSettings.tWindow(end,2) ];
    elseif strcmp(version, 'simulation')
        timeWindow  = [ TargetSettings.tIni, TargetSettings.tWindow(1,1) ];
    elseif strcmp(version, 'declustering')
        timeWindow  = [datetime(TargetSettings.startYear,1,1), datetime(9999,1,1)];
    else
        error('Unknown version')
    end
    isOutsideTimeWindow         	= Catalog.date < timeWindow(1) ...
                                    | Catalog.date > timeWindow(2);
    Catalog(isOutsideTimeWindow, :) = []; % & Catalog.mag < 5.5
    
    %% Filter by spatial window
    % Cut out events outside of complementary polygon
    isInPoly                = inpolygon( Catalog.lon, Catalog.lat, ...
                                         TargetSettings.polygonComple(:,1), TargetSettings.polygonComple(:,2) );
    Catalog(~isInPoly, :)   = [];
    
    %% Filter by magnitude & depth
    % Cut out events below/above threshold value
    Catalog( Catalog.mag < TargetSettings.Mc, : )   = [];
    Catalog( Catalog.depth > TargetSettings.maxDepth, :) = [];
    
    %% For simulation only: Include event ID that is forced to sample
    if strcmp(version, 'simulation')
        if ~isempty(TargetSettings.eventID_forcedSample)
            if ~ismember( TargetSettings.eventID_forcedSample, Catalog.id )
                Catalog = [Catalog; origCatalog( origCatalog.id==TargetSettings.eventID_forcedSample, : )];
            end
        end
    end
                               
end

%% Old code pieces
% complRangeX         = [min(polygon(:,1))-1, max(polygon(:,1))+1];
% complRangeY         = [min(polygon(:,2))-1, max(polygon(:,2))+1];
% polygon             = [ complRangeX(1), complRangeY(1); complRangeX(2), complRangeY(1); ...
%                         complRangeX(2), complRangeY(2); complRangeX(1), complRangeY(2) ];

