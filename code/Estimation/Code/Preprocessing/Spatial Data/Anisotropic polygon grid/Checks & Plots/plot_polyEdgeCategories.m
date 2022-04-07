function plot_anisoSpatPrecomp( Catalog, iEvent, origPolygon, rotatedPolygon )

% Convert strike to mathematical rotation angle
alpha       = 90-Catalog.strike(iEvent);
% Compute rupture start and end point
rupStart    = [Catalog.x(iEvent), Catalog.y(iEvent)] + Catalog.rupExtent(iEvent)/2 * [ -cosd(alpha), -sind(alpha) ];
rupEnd      = [Catalog.x(iEvent), Catalog.y(iEvent)] + Catalog.rupExtent(iEvent)/2 * [ cosd(alpha), sind(alpha) ];
                                
%% Event location and original polygon
figure;
subplot(1,3,1)
plot(origPolygon(1,:),origPolygon(2,:), 'k-')
hold on
plot(Catalog.x(iEvent), Catalog.y(iEvent), 'r*')
plot([rupStart(1), rupEnd(1)], [rupStart(2), rupEnd(2)], 'r-') 
axis([-10 10 -10 10])
axis square
title({'Event location and polygon before rotation', ['(strike = ', num2str(Catalog.strike(iEvent)), ')']})

%% Rotated event location and rotated polygon with colors indicating polygon edge category
subplot(1,3,2)
for i=1:size(rotatedPolygon,2)-1
    if flag(i)==-1
        color='r-'; 
    elseif flag(i)==1
        color='g-'; 
    elseif flag(i)==2
        color='b-';
    elseif flag(i)==3
        color='m';
    else
        color='k-';
    end
    plot(rotatedPolygon(1,[i:i+1]),rotatedPolygon(2,[i:i+1]), color, 'LineWidth', 2)
    hold on
end
plot(locR(1), locR(2), 'b*')
plot(locR(1)+[-rupLength/2 rupLength/2], repelem(locR(2),1,2), 'k-') 
axis([-10 10 -10 10])
axis square
title(['Categorization of polygon edges after rotation by ', num2str((strike-90)), ' degrees'])
