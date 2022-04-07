function [polygon, distances] = create_circularPolygon( radius_km, ...
                                                        centerLon, ...
                                                        centerLat )
% Call: [polygon, distances] = create_circularPolygon( 75, Catalog.lon(220235), Catalog.lat(220235) );
% [polygonTarget, distances] = create_circularPolygon( 75, Catalog.lon(220235), Catalog.lat(220235) );
% max(abs(distances-75))
% [polygonComple, distances] = create_circularPolygon( 100, Catalog.lon(220235), Catalog.lat(220235) );
% max(abs(distances-100))
% plot_temporalEventOccurrences( Catalog, polygonTarget, polygonComple )
% polygonTarget = create_circularPolygon( 75, Catalog.lon(648843), Catalog.lat(648843) );

    gridSize    = 101; % must be odd number
    if mod(gridSize,2)==0
        error('gridSize must be odd number')
    end
    distOneLat  = getDistance(0,0,0,1,'km');
    gridLat     = linspace(centerLat-radius_km/distOneLat, centerLat+radius_km/distOneLat, gridSize)';
    gridLon     = centerLon * ones(length(gridLat),1);
    lastDistance= 0;
    newLon      = centerLon;
    iterLon     = (gridLat(2)-gridLat(1))/100;
    
    % Find longitude corresponding to latitude
    for i=2:length(gridLat)-1
        
        currentDistance = getDistance( newLon, gridLat(i), centerLon, centerLat, 'km' );
        
        if i <= ceil(gridSize/2)
            
            while currentDistance < radius_km                    
                lastLon         = newLon;
                newLon          = lastLon + iterLon;
    %             iLonDiff        = iLonDiff + iterLonDiff;
                lastDistance    = currentDistance;
                currentDistance = getDistance( newLon, gridLat(i), centerLon, centerLat, 'km' );
            end
            
        else
            
            while currentDistance > radius_km                    
                lastLon         = newLon;
                newLon          = lastLon - iterLon;
    %             iLonDiff        = iLonDiff + iterLonDiff;
                lastDistance    = currentDistance;
                currentDistance = getDistance( newLon, gridLat(i), centerLon, centerLat, 'km' );
            end
            
        end
        
        factor      = (radius_km - lastDistance) / (currentDistance - lastDistance);
        gridLon(i)  = lastLon + factor * (newLon-lastLon);
        
    end
    
    gridLat2 = flipud( gridLat(1:end-1) );
    gridLon2 = flipud( centerLon - (gridLon(1:end-1) - centerLon) );
    
    polygon = [gridLon, gridLat; gridLon2, gridLat2];
    distances = getDistance( polygon(:,1), polygon(:,2), centerLon, centerLat, 'km' );
    
    % Make anti-clockwise
    
end