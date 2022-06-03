function [ gridX, gridY ] = create_gridXY4polygon( polygon )

    % x-y ranges
    xRange = [min(polygon(:,1)), max(polygon(:,1))];
    yRange = [min(polygon(:,2)), max(polygon(:,2))];
  
    % Grid dimensions
    rv = diff(yRange)/diff(xRange);
    dim_xy = round( 128 * ( [1, rv]*(rv > 1) + [1/rv, 1]*(rv <= 1) ) );

    % Meshgrid
    gridX                       = linspace(xRange(1), xRange(2), dim_xy(1));
    gridY                       = linspace(yRange(1), yRange(2), dim_xy(2));
%     dX                          = abs(gridX(2)-gridX(1));
%     dY                          = abs(gridY(2)-gridY(1));
    [ gridY_mesh, gridX_mesh ]  = meshgrid( gridY, gridX );
    gridX                       = gridX_mesh(:)';
    gridY                       = gridY_mesh(:)';
    
    % Cut out mesh grid points outside of polygon
    isInPolygon                 = inpolygon( gridX, gridY, polygon(:,1), polygon(:,2) );
    gridX( ~isInPolygon )       = [];
    gridY( ~isInPolygon )       = [];
    
end