function spatVolatility = evaluate_spatialClustering( Catalog, isTarget, catalogID, polygon, polyArea )

    %     if catalogID <= 0
%         [ gridX, gridY ] = create_gridXY4polygon( polygon );
%         dist            = getDistance(Catalog.x(isTarget), Catalog.y(isTarget), gridX, gridY, spaceUnit);
%         [~, gridIdx]    = min(dist,[], 2);
%         hist_spat       = histcounts(gridIdx, 0.5:1:length(gridX)+0.5)';
%         spatVolatility  = std(hist_spat)/mean(hist_spat);
%     else
%         spatVolatility  = 0;
%     end
    if catalogID <= 0
        spatVolatility  = {estimate_ripleyK( Catalog.x(isTarget), Catalog.y(isTarget), polyArea )};
    else
        spatVolatility = {0};
    end
    
end