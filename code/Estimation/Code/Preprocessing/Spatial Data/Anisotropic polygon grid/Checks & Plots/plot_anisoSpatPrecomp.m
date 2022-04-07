function plot_anisoSpatPrecomp( AnisoGrid, IsoLeftGrid, IsoRightGrid, Catalog, iEvent, Rotated, tMinMaxDist  )

    % Convert strike to mathematical rotation angle
    alpha       = 90-Catalog.strike(iEvent);
    % Compute rupture start and end point
    rupStart    = [Catalog.x(iEvent), Catalog.y(iEvent)] + Catalog.rupExtent(iEvent)/2 * [ -cosd(alpha), -sind(alpha) ];
    rupEnd      = [Catalog.x(iEvent), Catalog.y(iEvent)] + Catalog.rupExtent(iEvent)/2 * [ cosd(alpha), sind(alpha) ];
    
    % If structure is empty, fill it with empty variables
    if isempty(fieldnames(AnisoGrid))
        AnisoGrid.x             = [];
        AnisoGrid.y             = [];
        AnisoGrid.r             = [];
        AnisoGrid.polyStart     = [];
        AnisoGrid.polyEnd       = [];
        AnisoGrid.polyCategory  = [];
        AnisoGrid.segmWeight    = [];
        AnisoGrid.segmFactor    = [];
    end
    if isempty(fieldnames(IsoLeftGrid))
        IsoLeftGrid.x             = [];
        IsoLeftGrid.y             = [];
        IsoLeftGrid.r             = [];
        IsoLeftGrid.polyStart     = [];
        IsoLeftGrid.polyEnd       = [];
        IsoLeftGrid.polyCategory  = [];
        IsoLeftGrid.segmWeight    = [];
        IsoLeftGrid.segmFactor    = [];
    end
    if isempty(fieldnames(IsoRightGrid))
        IsoRightGrid.x             = [];
        IsoRightGrid.y             = [];
        IsoRightGrid.r             = [];
        IsoRightGrid.polyStart     = [];
        IsoRightGrid.polyEnd       = [];
        IsoRightGrid.polyCategory  = [];
        IsoRightGrid.segmWeight    = [];
        IsoRightGrid.segmFactor    = [];
    end

    %% Event location and original polygon
    figure;
    subplot(1,3,1)
    plot(Rotated.origPolygon(1,:),Rotated.origPolygon(2,:), 'k-')
    hold on
    plot(Catalog.x(iEvent), Catalog.y(iEvent), 'r*')
    plot([rupStart(1), rupEnd(1)], [rupStart(2), rupEnd(2)], 'r-') 
    axis([-10 17 -10 17])
    axis square
    title({'Event location and polygon before rotation', ['(strike = ', num2str(Catalog.strike(iEvent)), ')']})

    %% Rotated polygon with colors indicating postive/negative contribution to integral
    subplot(1,3,2)
    l1 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', 'r',           'LineWidth', 2);
    hold on
    l2 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', 'g',           'LineWidth', 2);
    l3 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', [0.6 0 0.2],   'LineWidth', 2);
    l4 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', [1 0.7 0.1],   'LineWidth', 2);
    l5 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', [0.75 0 0.75], 'LineWidth', 2);
    l6 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', [0.3 0.7 0.9], 'LineWidth', 2);
    
    % Plot anisotropic parts
    for i=1:size(AnisoGrid.polyStart,1)
        linestyle   = '-';
        if AnisoGrid.polyCategory(i) < 1
            color       = 'r';            
        else
            color       = 'g';
        end
        lineX = [AnisoGrid.polyStart(i,1), AnisoGrid.polyEnd(i,1)];
        lineY = [AnisoGrid.polyStart(i,2), AnisoGrid.polyEnd(i,2)];
        plot(lineX, lineY, 'LineStyle', linestyle, 'Color', color, 'LineWidth', 2)
        hold on
    end
    % Plot isotropic left parts
    for i=1:size(IsoLeftGrid.polyStart,1)
        linestyle   = '-';
        if IsoLeftGrid.polyCategory(i) < 1
            color       = [0.6 0 0.2];            
        else
            color       = [1 0.7 0.1];
        end
        lineX = [IsoLeftGrid.polyStart(i,1), IsoLeftGrid.polyEnd(i,1)];
        lineY = [IsoLeftGrid.polyStart(i,2), IsoLeftGrid.polyEnd(i,2)];
        plot(lineX, lineY, 'LineStyle', linestyle, 'Color', color, 'LineWidth', 2)
        hold on
    end
    % Plot isotropic right parts
    for i=1:size(IsoRightGrid.polyStart,1)
        linestyle   = '-';
        if IsoRightGrid.polyCategory(i) < 1
            color       = [0.75 0 0.75];            
        else
            color       = [0.3 0.7 0.9];
        end
        lineX = [IsoRightGrid.polyStart(i,1), IsoRightGrid.polyEnd(i,1)];
        lineY = [IsoRightGrid.polyStart(i,2), IsoRightGrid.polyEnd(i,2)];
        plot(lineX, lineY, 'LineStyle', linestyle, 'Color', color, 'LineWidth', 2)
        hold on
    end
    plot(Rotated.eventLoc(1), Rotated.eventLoc(2), 'r*')
    plot([Rotated.rupXleft, Rotated.rupXright], [Rotated.rupY, Rotated.rupY], 'k-', 'LineWidth', 2) 
    l7 = plot([Rotated.rupXleft, Rotated.rupXleft], [-10, 15], 'k', 'LineStyle', ':', 'LineWidth', 0.1);
    plot([Rotated.rupXright, Rotated.rupXright], [-10, 15], 'k', 'LineStyle', ':', 'LineWidth', 0.1)
    axis([-10 17 -10 17])
    axis square
    title({'Rotated event location and polygon', ...
           'Colors indicate neg./pos. contribution to integral' ...
         })
    xlabel({['Sum of weights: Aniso= ', num2str(sum(AnisoGrid.segmWeight)), '; IsoLeft= ', num2str(sum(IsoLeftGrid.segmWeight)), '; IsoRight= ', num2str(sum(IsoRightGrid.segmWeight))], ...
            ['Sum of factors: Aniso= ', num2str(sum(AnisoGrid.segmFactor)), '; IsoLeft= ', num2str(sum(IsoLeftGrid.segmFactor)), '; IsoRight= ', num2str(sum(IsoRightGrid.segmFactor))] ...
          })
    legend([l1 l2 l3 l4 l5 l6 l7], {'- aniso', '+ aniso', '- isoLeft', '+ isoLeft', '- isoRight', '+ isoRight', 'aniso area'})
    
    %% Rotated polygon with colors indicating distance to rupture line
    subplot(1,3,3)
    x = [ AnisoGrid.x; IsoLeftGrid.x; IsoRightGrid.x ];
    y = [ AnisoGrid.y; IsoLeftGrid.y; IsoRightGrid.y ];
    r = [ AnisoGrid.r; IsoLeftGrid.r; IsoRightGrid.r ];
    D       = 0.01;
    gamma   = 1;
    q       = 2;
    minDist = tMinMaxDist.r_min(iEvent);
    maxDist = tMinMaxDist.r_max(iEvent);
    F_min   = 1 - (1 + (2*Catalog.rupExtent(iEvent).*minDist + pi*minDist^2)./(D*exp(gamma*Catalog.mag(iEvent))))^(1-q);
    warning('Magnitude needs to be reduced by magnitude threshold!!')
    F_max   = 1 - (1 + (2*Catalog.rupExtent(iEvent).*maxDist + pi*maxDist^2)./(D*exp(gamma*Catalog.mag(iEvent))))^(1-q);
    
    scatter(x, y, 3, r, 'filled')
    colorbar;
    hold on
    plot(Rotated.eventLoc(1), Rotated.eventLoc(2), 'r*')
    plot([Rotated.rupXleft, Rotated.rupXright], [Rotated.rupY, Rotated.rupY], 'r-', 'LineWidth', 2) 
    axis([-10 17 -10 17])
    axis square
    title({'Distances of polygon edge grid points to', 'event rupture (in units of latitude degrees)'})
    % Event in polygon
    if abs(Catalog.flag(iEvent)) >= 0.5
        xlabel({['Min. distance: ', num2str(minDist), ' => integral over this distance: ', num2str(F_min)], ...
                ['Max. distance: ', num2str(maxDist), ' => integral over this distance: ', num2str(F_max)], ...
                 '(initial guesses: D=0.01, gamma=1, q=2;', ...
                 '1 - (1 + (2*rupExtent*dist + pi*dist^2)./(D*exp(gamma*Magn)))^{(1-q)}' ...
               })
    else
        xlabel({['Min. distance: ', num2str(minDist)], ...
                ['Max. distance: ', num2str(maxDist)], ...
                 '(integral cannot be approximated since event lies outside of polygon)' ...
               })
    end
      
end
