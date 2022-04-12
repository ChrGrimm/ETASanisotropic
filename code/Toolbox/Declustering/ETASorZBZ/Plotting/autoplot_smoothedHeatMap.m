function autoplot_smoothedHeatMap( data )

    xlims = [floor(min(data(:,1))), ceil(max(data(:,1)))];
    ylims = [floor(min(data(:,2))), ceil(max(data(:,2)))];
    
    xIncr = xlims(1):0.02:xlims(2);
    yIncr = ylims(1):0.02:ylims(2);
    
    [X,Y] = meshgrid(xIncr, yIncr);
    nGrid = length(X(:));
    nOcc  = zeros(nGrid,1);
    
    for i=1:nGrid
        nOcc(i) = sum( getDistance(X(i), Y(i), data(:,1), data(:,2), 'degree') <= 0.1 );
    end
    
    scatter(X(:), Y(:), 5, nOcc, 'filled')
    hold on
    for i=1:6
        plot(xlims, -(xlims+i), 'r--')
    end
    axis([xlims, ylims])
    colorbar

end