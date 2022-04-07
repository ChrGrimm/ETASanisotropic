function plot_isotropicGrid( grid_x, ...
                             grid_y, ...
                             grid_r, ...
                             locX, ...
                             locY, ...
                             polyStart, ...
                             polyEnd, ...
                             segmWeight, ...
                             segmFactor, ...
                             signPolyEdge, ...
                             version )
         
    rangeX_left     = min([polyStart(:,1); polyEnd(:,1)]) -4;
    rangeX_right    = max([polyStart(:,1); polyEnd(:,1)]) +10;
    rangeY_bottom   = min([polyStart(:,2); polyEnd(:,2)]) -4;
    rangeY_top      = max([polyStart(:,2); polyEnd(:,2)]) +10;    
                         
    if strcmp(version, 'isotropic')
        
        close all
        figure;
        
        %% Rotated polygon with colors indicating postive/negative contribution to integral
        subplot(1,2,1)
        % Dummy plots for legend entries
        l1 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', 'r',           'LineWidth', 2);
        hold on
        l2 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', 'g',           'LineWidth', 2);

        % Plot polygon parts
        if ~isempty(polyStart)
            for i=1:size(polyStart,1)
                linestyle   = '-';
                if signPolyEdge(i) == -1
                    color       = 'r';            
                elseif signPolyEdge(i) == 1
                    color       = 'g';
                else
                    error('signPolyEdge must be other -1 or 1!')
                end
                lineX = [polyStart(i,1), polyEnd(i,1)];
                lineY = [polyStart(i,2), polyEnd(i,2)];
                plot(lineX, lineY, 'LineStyle', linestyle, 'Color', color, 'LineWidth', 2)
                hold on
            end
        end
        plot(locX, locY, 'r*')
        axis([rangeX_left, rangeX_right, rangeY_bottom, rangeY_top])
        axis square
        title({'Event location and polygon', ...
               'Colors indicate neg./pos. contribution to integral' ...
             })
        xlabel({['Sum of weights = ', num2str(sum(segmWeight))], ...
                ['Sum of factors = ', num2str(sum(segmFactor))] ...
              })
        legend([l1 l2], {'negative contribution', 'positive contribution'})

        %% Rotated polygon with colors indicating distance to rupture line
        subplot(1,2,2)

        scatter(grid_x, grid_y, 3, grid_r, 'filled')
        colorbar;
        hold on
        plot(locX, locY, 'r*') 
        axis([rangeX_left, rangeX_right, rangeY_bottom, rangeY_top])
        axis square
        title({'Distances of polygon edge grid points to', 'event location (in units of latitude degrees)'})
        
    elseif strcmp(version, 'isotropic left') || strcmp(version, 'isotropic right')
        
        subplot(1,3,2)
        hold on
        % Plot isotropic left parts
        for i=1:size(polyStart,1)
            linestyle   = '-';
            if signPolyEdge(i) == -1
                if strcmp(version, 'isotropic left')
                    color = [0.6 0 0.2];
                else
                    color = [0.75 0 0.75];
                end
            elseif signPolyEdge(i) == 1
                if strcmp(version, 'isotropic right')
                    color = [1 0.7 0.1];
                else
                    color = [0.3 0.7 0.9];
                end
            else
                error('Sign must be either -1 or 1!')
            end
            lineX = [polyStart(i,1), polyEnd(i,1)];
            lineY = [polyStart(i,2), polyEnd(i,2)];
            plot(lineX, lineY, 'LineStyle', linestyle, 'Color', color, 'LineWidth', 2)
            hold on
        end
        
        %%
        subplot(1,3,3)
        scatter(grid_x, grid_y, 3, grid_r, 'filled')
        
    end
                           
end