function [ polyRotated, ...
           polyStart, ...
           polyEnd, ...
           polyEdgeCategory, ...
           locRotated, ...
           rupLength, ...
           rupXleft, ...
           rupXright, ...
           rupY, ...
           rupFlag ] = categorize_polyEdges( Catalog, polygon, iEvent, makePlots )
    
%     % Initialize inputs
%     strike = 45;
%     polygon = [2,2; 0.5,1; -0.5,-0.5; -1.5,1; -0.8,-2; -0.8, -1.5; 2, 1; 2,2]';
%     loc = [0, 0]';
    
    %% Extract event information
    loc         = [ Catalog.x(iEvent); Catalog.y(iEvent) ];
    rupLength   = Catalog.rupExtent(iEvent);
    strike      = Catalog.strike(iEvent);
    
    %% Compute rotated location and polygon
    RotMatrix   = compute_rotationMatrix( strike, 'horizontal' );
    locRotated  = RotMatrix * loc;
    polyRotated = RotMatrix * polygon;

    %% Compute rupture coordinates
    rupXleft            = locRotated(1)-rupLength/2;
    rupXright           = locRotated(1)+rupLength/2;
    rupY                = locRotated(2);
    % Rupture flag: 
    % 1     if rupture end point lies within polygon
    % 0.5   if rupture end point lies on polygon edge
    % 0     if rupture end point lies outside of polygon
    [inPoly, onPoly]    = inpolygon([rupXleft, rupXright], [rupY, rupY], polyRotated(1,:), polyRotated(2,:));
    rupFlag             = inPoly - 0.5 * onPoly;

    %% Split and categorize polygon segments
    polyEdgeCategory = [];
    i = 1;
    while i < size(polyRotated,2)

       x1 = polyRotated(1,i);
       x2 = polyRotated(1,i+1);
       y1 = polyRotated(2,i);
       y2 = polyRotated(2,i+1);
       % polyEdgeCategory = -1: parallel to rupture line, negative orientation
       % polyEdgeCategory = 0: vertcial to or lying on rupture line
       % polyEdgeCategory = 1: parallel to rupture line, positive orientation      
       % polyEdgeCategory = 2: circle left (is set to -2 in compute_isotropicGrid, if negative orientation)
       % polyEdgeCategory = 3: circle right (is set to -3 in compute_isotropicGrid, if negative orientation)
       % polygon edge outside of rupture box
       if max(x1,x2) <= rupXleft
           polyEdgeCategory = [polyEdgeCategory, 2];
           i = i + 1;
       elseif min(x1,x2) >= rupXright
           polyEdgeCategory = [polyEdgeCategory, 3];
           i = i + 1;
       % polygone edge inside of rupture box
       elseif min(x1,x2) >=  rupXleft && max(x1,x2) <= rupXright
           % polygon edge vertical to rupture line or lying on rupture line
           if (x1 == x2) || (y1 == rupY && y2 == rupY)
               polyEdgeCategory = [polyEdgeCategory, 0];
               i = i + 1;
           % polygon edge oriented towards rupture line
           elseif (min(y1,y2)>=rupY && x2<x1) || (max(y1,y2)<=rupY && x2>x1)
               polyEdgeCategory = [polyEdgeCategory, 1];
               i = i + 1;
           % polygon edge oriented away from rupture line
           elseif (min(y1,y2)>=rupY && x2>x1) || (max(y1,y2)<=rupY && x2<x1)
               polyEdgeCategory = [polyEdgeCategory, -1];
               i = i + 1;
           % polygone crosses rupture line
           % elseif (y1>rupY && y2<rupY) || (y1<rupY && y2>rupY)
           else
               xnew = x1+(x2-x1)*(rupY-y1)/(y2-y1);
               polyRotated = [polyRotated(:,1:i), [xnew; rupY], polyRotated(:,(i+1):end)];
           end
       % polygone edge crossing vertical boundaries of rupture box
       else
           % crosses rupture box on right hand side
           if (x1>rupXright && x2<rupXright) || (x1<rupXright && x2>rupXright)
               ynew = y1+(y2-y1)*(rupXright-x1)/(x2-x1);
               polyRotated = [polyRotated(:,1:i), [rupXright; ynew], polyRotated(:,(i+1):end)];
           % crosses rupture box on right hand side
           % elseif (x1>rupXleft && x2<rupXleft) || (x1<rupXleft && x2>rupXleft)
           else
               ynew = y1+(y2-y1)*(rupXleft-x1)/(x2-x1);
               polyRotated = [polyRotated(:,1:i), [rupXleft; ynew], polyRotated(:,(i+1):end)];
           end   
       end

    end

%     %% Store data in struct
%     Rotated = struct('polygon', polyRotated, 'polyEdgeCategory', polyEdgeCategory, 'rupLength', rupLength, ...
%                      'rupXleft', rupXleft, 'rupXright', rupXright, 'rupY', rupY, 'eventLoc', locRotated);
                 
    % Start and end points of rotated polygon edges
    polyStart   = polyRotated(:,1:end-1);
    polyEnd     = polyRotated(:,2:end);

    if makePlots
        
        close all
        
        %% Event location and original polygon
        figure;
        subplot(1,3,1)
        plot(polygon(1,:),polygon(2,:), 'k-')
        hold on
        plot(Catalog.x(iEvent), Catalog.y(iEvent), 'r*')
        origRupLeft     = loc - rupLength/2 * [ cosd(90-strike); sind(90-strike) ];
        origRupRight    = loc + rupLength/2 * [ cosd(90-strike); sind(90-strike) ];
        plot([origRupLeft(1), origRupRight(1)], [origRupLeft(2), origRupRight(2)], 'r-') 
        axis([-10 17 -10 17])
        axis square
        title({'Event location and polygon before rotation', ['(strike = ', num2str(Catalog.strike(iEvent)), ')']})
        
        subplot(1,3,2)
        l1 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', 'r',           'LineWidth', 2);
        hold on
        l2 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', 'g',           'LineWidth', 2);
        l3 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', [0.6 0 0.2],   'LineWidth', 2);
        l4 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', [1 0.7 0.1],   'LineWidth', 2);
        l5 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', [0.75 0 0.75], 'LineWidth', 2);
        l6 = plot(nan(2,1), nan(2,1), 'LineStyle', '-', 'Color', [0.3 0.7 0.9], 'LineWidth', 2);
        l7 = plot(nan(2,1), nan(2,1), 'k', 'LineStyle', ':', 'LineWidth', 0.1);
        
        plot(locRotated(1), locRotated(2), 'r*')
        plot([rupXleft, rupXright], [rupY, rupY], 'r-', 'LineWidth', 2)
        
        legend([l1 l2 l3 l4 l5 l6 l7], ...
               {'- aniso', '+ aniso', '- isoLeft', '+ isoLeft', '- isoRight', '+ isoRight', 'aniso area'}, ...
               'AutoUpdate','off')
        axis([-10 17 -10 17])
        axis square
        title({'Rotated event location and polygon', ...
               'Colors indicate neg./pos. contribution to integral' ...
             })
         
        %%
        subplot(1,3,3)
        plot(locRotated(1), locRotated(2), 'r*')
        hold on
        plot([rupXleft, rupXright], [rupY, rupY], 'r-', 'LineWidth', 2) 
        axis([-10 17 -10 17])
        axis square
        title({'Distances of polygon edge grid points to', 'event rupture (in units of latitude degrees)'})
         
    end
    
end