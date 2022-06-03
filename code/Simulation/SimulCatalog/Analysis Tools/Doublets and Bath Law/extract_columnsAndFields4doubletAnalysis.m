function [ xCoord, yCoord, t, mag, date, flag, ...
           eventID, triggerID, clusterID, ...
           timeWindow_days, spatialWindow, ...
           magnWindow, magnIntervals, nMagnitudes ] = extract_columnsAndFields4doubletAnalysis( Catalog, ...
                                                                                                DoubletCriteria, ...
                                                                                                isSynthetic )
                                                                                            
    %% Extract columns from table Catalog
    xCoord          = Catalog.x;
    yCoord          = Catalog.y;
    t               = Catalog.t;
    mag             = round(Catalog.mag,1);
    flag            = Catalog.flag;
    eventID         = Catalog.id;
    triggerID       = Catalog.triggerID;
    clusterID       = Catalog.clusterID;
    if isSynthetic
        date        = NaN;    
    else
        date        = Catalog.date;
    end
    
    %% Extract fields from struct DoubletCriteria
    % Time window
    timeWindow_days = DoubletCriteria.timeWindow_days;       

    % Evaluate event-specific spatial window
    spatialWindow    = max(DoubletCriteria.factorSpatExtent * Catalog.rupL, DoubletCriteria.minSpatExtent);
    
    % Magnitude window and intervals (round magnitudes)
    magnWindow      = round(DoubletCriteria.magnWindow,1);
    magnIntervals   = round(DoubletCriteria.magnIntervals,1);
    nMagnitudes     = length(magnIntervals);
    
                                                                                            
end