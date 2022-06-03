function [ spatialDistrib, ...
           nEventsPerGridCell ] = evaluate_spatialEventDistrib( xEvents, yEvents, xGrid, yGrid )

    idxNearestGridCell  = knnsearch( [xGrid,yGrid], [xEvents, yEvents] );
    nEventsPerGridCell  = hist( idxNearestGridCell, 1:length(xGrid) );
    nEventsTotal        = length(xEvents);
    spatialDistrib      = nEventsPerGridCell / nEventsTotal;

end