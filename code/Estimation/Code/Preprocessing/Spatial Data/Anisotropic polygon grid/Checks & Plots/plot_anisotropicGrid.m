function plot_anisotropicGrid( grid_x, ...
                               grid_y, ...
                               grid_r, ...
                               rupXleft, ...
                               rupXright, ...
                               rupY, ...
                               polyStart, ...
                               polyEnd, ...
                               polyEdgeCategory )
    
    subplot(1,3,2)
    hold on
            
    plot([rupXleft, rupXleft], [-10, 15], 'k', 'LineStyle', ':', 'LineWidth', 0.1);
    plot([rupXright, rupXright], [-10, 15], 'k', 'LineStyle', ':', 'LineWidth', 0.1)

    % Plot anisotropic parts
    for i=1:size(polyStart,1)
        linestyle   = '-';
        if polyEdgeCategory(i) < 1
            color       = 'r';            
        else
            color       = 'g';
        end
        lineX = [polyStart(i,1), polyEnd(i,1)];
        lineY = [polyStart(i,2), polyEnd(i,2)];
        plot(lineX, lineY, 'LineStyle', linestyle, 'Color', color, 'LineWidth', 2)
    end
    
    %%
    subplot(1,3,3)
    hold on
    scatter(grid_x, grid_y, 3, grid_r, 'filled')
    colorbar;
                           
end