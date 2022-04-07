function [ x_central, ...
           y_central, ...
           polygon_central, ...
           polygonCenter, ...
           polyArea ] = convert_lonlat2xy( lon, ...
                                             lat, ...
                                             polygon )
                                                         
    % Calculate centroid of polygon region
    polyin                      = polyshape( polygon(:,1), polygon(:,2) );
    [polyCenterX, polyCenterY]  = centroid( polyin );
    % Latitude factor
    latitudeFactor              = cos(polyCenterY * pi / 180);
    % Centralize events' geo-coordinates and apply latitude factor
    x_central                   = (lon - polyCenterX) * latitudeFactor;
    y_central                   = lat - polyCenterY;
    % Centralize polygons' geo-coordinates and apply latitude factor
    polygon_central(:,1)        = (polygon(:,1) - polyCenterX) * latitudeFactor;
    polygon_central(:,2)        = polygon(:,2) - polyCenterY;
    % Save polygon center coordinates in struct TargetWindow
    polygonCenter               = [polyCenterX, polyCenterY];
    
    %% Polygon area
    polyArea                    = polyarea( polygon_central(:,1), polygon_central(:,2) );

end