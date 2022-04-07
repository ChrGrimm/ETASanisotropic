function plot_anisoEventDistr( Catalog, iEvent, f_inn, spatParamETAS )

    D = spatParamETAS(1);
    gamma = spatParamETAS(2);
    q = spatParamETAS(3);

    polygon = [ -6.076 6.076 6.076 -6.076 -6.076; -8 -8 8 8 -8 ];
    % Convert strike to mathematical rotation angle
    alpha       = 90-Catalog.strike(iEvent);
    % Compute rupture start and end point
    rupStart    = [Catalog.x(iEvent), Catalog.y(iEvent)] + Catalog.rupExtent(iEvent)/2 * [ -cosd(alpha), -sind(alpha) ];
    rupEnd      = [Catalog.x(iEvent), Catalog.y(iEvent)] + Catalog.rupExtent(iEvent)/2 * [ cosd(alpha), sind(alpha) ];
    
    %% Event location and original polygon
    figure;
    subplot(1,2,1)
    plot(polygon(1,:),polygon(2,:), 'k-')
    hold on
    plot(Catalog.x(iEvent), Catalog.y(iEvent), 'r*')
    plot([rupStart(1), rupEnd(1)], [rupStart(2), rupEnd(2)], 'r-') 
    axis([-10 17 -10 17])
    axis square
    title({'Event location and polygon before rotation', ['(strike = ', num2str(Catalog.strike(iEvent)), ')']})

    %% Event spatial distribution
    subplot(1,2,2)
    [X,Y] = meshgrid(-2*Catalog.rupExtent(iEvent):0.01:2*Catalog.rupExtent(iEvent));
    r = get_dist2line([-Catalog.rupExtent(iEvent)/2, 0], [ Catalog.rupExtent(iEvent)/2, 0], [X(:),Y(:)]', 0)';
    f = (q-1)/(D*exp(gamma*Catalog.mag(iEvent))) * f_inn(r, r.^2, Catalog.rupExtent(iEvent), Catalog.mag(iEvent)).^(-q);
    warning('Magnitude needs to be reduced by magnitude threshold!!')
    f = reshape(f, size(X));
    h = surf(X, Y, f);
%         set(h,'LineStyle','none')
    title(['PDF spatial distribution for Event ', num2str(iEvent), ' with Mw = ', num2str(Catalog.mag(iEvent))])
end