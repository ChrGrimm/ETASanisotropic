[magSorted,idx] = sort(Catalog.mag, 'descend');
CatalogStrongest = Catalog(idx(1:10),:);
CatalogStrongest.origStr1 = tempCatalog.str1(idx(1:10));
CatalogStrongest.origStr2 = tempCatalog.str2(idx(1:10));
CatalogStrongest = sortrows(CatalogStrongest, 'date');

figure
for i=1:10
    subplot(2,5,i)
    iStrike = CatalogStrongest.strike(i);
    t       = Catalog.t;
    iT      = CatalogStrongest.t(i);
    x       = Catalog.x;
    iX      = CatalogStrongest.x(i);
    y       = Catalog.y;
    iY      = CatalogStrongest.y(i);
    iRupL   = CatalogStrongest.rupL(i);
    iEpiPos = CatalogStrongest.epiPos(i);
    RotMatrix   = compute_rotationMatrix( iStrike, 'vertical' );
    isAftersh   = t>iT & t<=iT+2 & getDistance(x, y, iX, iY, Method.spaceUnit) <= 2.5*iRupL;
    locsAftersh = [x(isAftersh)'; y(isAftersh)'];
    locsAftershRotated = RotMatrix * locsAftersh;
    locsMainshRotated = RotMatrix * [iX; iY];

    % plot_shapeMap
    plot(locsMainshRotated(1), locsMainshRotated(2), 'r*', 'LineWidth', 2)
    hold on
    scatter(locsAftershRotated(1,:), locsAftershRotated(2,:), 5, 'blue', 'filled')
    plot(locsMainshRotated(1)+[0,0], ...
          locsMainshRotated(2)+iRupL*[-iEpiPos, 1-iEpiPos], ...
          'r-', 'LineWidth', 2)
    origStrRot1 = CatalogStrongest.origStr1(i)+iStrike;
    origStrRot2 = CatalogStrongest.origStr2(i)+iStrike;
    plot(locsMainshRotated(1)+iRupL/2*cos(origStrRot1)*[-1,1], ...
          locsMainshRotated(2)+iRupL/2*sin(origStrRot1)*[-1,1], ..., ...
          'g-', 'LineWidth', 1.5)
    plot(locsMainshRotated(1)+iRupL/2*cos(origStrRot2)*[-1,1], ...
          locsMainshRotated(2)+iRupL/2*sin(origStrRot2)*[-1,1], ...
          'g-', 'LineWidth', 1.5)
      
    title([datestr(CatalogStrongest.date(i)), ', Mw=', num2str(CatalogStrongest.mag(i))])
    
end
