function plot_isoSpatPrecomp( IsoGrid, Catalog, iEvent )

    %% Rotated polygon with colors indicating postive/negative contribution to integral
    subplot(1,2,1)
    % Dummy plots for legend entries
    l1 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', 'r',           'LineWidth', 2);
    hold on
    l2 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', 'g',           'LineWidth', 2);
    
    % Plot polygon parts
    if ~isempty(fieldnames(IsoGrid))
        for i=1:size(IsoGrid.polyStart,1)
            linestyle   = '-';
            if IsoGrid.polyCategory(i) < 1
                color       = 'r';            
            else
                color       = 'g';
            end
            lineX = [IsoGrid.polyStart(i,1), IsoGrid.polyEnd(i,1)];
            lineY = [IsoGrid.polyStart(i,2), IsoGrid.polyEnd(i,2)];
            plot(lineX, lineY, 'LineStyle', linestyle, 'Color', color, 'LineWidth', 2)
            hold on
        end
    end
    plot(Catalog.x(iEvent), Catalog.y(iEvent), 'r*')
    axis([-10 17 -10 17])
    axis square
    title({'Event location and polygon', ...
           'Colors indicate neg./pos. contribution to integral' ...
         })
    xlabel({['Sum of weights = ', num2str(sum(IsoGrid.segmWeight))], ...
            ['Sum of factors = ', num2str(sum(IsoGrid.segmFactor))] ...
          })
    legend([l1 l2], {'negative contribution', 'positive contribution'})
    
    %% Rotated polygon with colors indicating distance to rupture line
    subplot(1,2,2)
    
    scatter(IsoGrid.x, IsoGrid.y, 3, IsoGrid.r, 'filled')
    colorbar;
    hold on
    plot(Catalog.x(iEvent), Catalog.y(iEvent), 'r*') 
    axis([-10 17 -10 17])
    axis square
    title({'Distances of polygon edge grid points to', 'event location (in units of latitude degrees)'})
      
end
