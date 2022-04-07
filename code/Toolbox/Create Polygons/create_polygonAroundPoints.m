function polygon = create_polygonAroundPoints( x, y, shrinkFactor )

    k = boundary(x, y, shrinkFactor);
    polygon(:,1) = x(k);
    polygon(:,2) = y(k);
    
    figure
    scatter(x, y, 2, 'b', 'filled')
    hold on
    plot(polygon(:,1), polygon(:,2), 'r-')

end