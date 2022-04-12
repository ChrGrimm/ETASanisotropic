function newCluster = identify_eventsInSpaceTimeWindow( WindowSettings, Catalog, iTriggEvent )

    %% Initialize output table
    newCluster = table;
    
    %% Apply time limit
    timeDiff        = days( Catalog.date - Catalog.date(iTriggEvent) );
    idxInTime       = find( timeDiff >= 0 & timeDiff <= WindowSettings.tWindow_days );

    %% Apply space limit
    if Catalog.lon(iTriggEvent) < -90 || Catalog.lon(iTriggEvent) > 90
        Catalog.lon = correct_longitudesForDateLine( Catalog.lon );
    end
    
    distance        = getDistance( Catalog.lon(idxInTime), Catalog.lat(idxInTime), ...
                                   Catalog.lon(iTriggEvent), Catalog.lat(iTriggEvent), 'km' );
    idxInTimeSpace  = idxInTime( distance <= Catalog.spatRestr(iTriggEvent) );
    
    %% Save results in table
    newCluster.evIDs = Catalog.id(idxInTimeSpace);

end