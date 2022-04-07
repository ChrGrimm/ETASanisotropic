if strcmp(Method.spaceUnit, 'degree')
    % Decentralize event locations
    lon = Catalog.x / cos(TargetSettings.polyCenter(2) * pi / 180) + TargetSettings.polyCenter(1);
    lat = Catalog.y + TargetSettings.polyCenter(2);
else
    lon = Catalog.x;
    lat = Catalog.y;
end

magnPlus = Catalog.mag(eventID);

%% Mesh grid and distances to event/ rupture line
[rupExtent, rupW] = estimate_rupSize(TargetSettings.M_c + magnPlus, Tectonics.type, Tectonics.faultingStyle, Method.spaceUnit);
k               = Method.factorSpatialRestr;
c               = Method.minSpatialRestr;
spatialExtent   = max(k * rupExtent, c);
limitMesh = 3*rupExtent; %max(limitMesh, 1.1*spatialExtent);
[X,Y]   = meshgrid(-limitMesh:limitMesh/100:limitMesh);
X       = X + lon(eventID);
Y       = Y + lat(eventID);

if startsWith(Method.spatialDesign, 'aniso')
    
    [rupStart, rupEnd] = compute_rupCoords( lon(eventID), ...
                                            lat(eventID), ...
                                            Catalog.rupExtent(eventID), ...
                                            Catalog.str(eventID), ...
                                            Catalog.epiPos(eventID) );
                                        
    r = get_dist2line(rupStart, rupEnd, [X(:),Y(:)]', Method.spaceUnit, 0)';
    
else
    
    r = sqrt( (X-lon(eventID)).^2 + (Y-lat(eventID)).^2 );
    
end

f = spatialFactor(spatialExtent,magnPlus,rupExtent,rupW,dip) .* f_inn(r,r.^2, magnPlus, rupExtent, rupW, dip).^(1-q) .* (r <= spatialExtent);

plot_shapeMap
scatter(X(:), Y(:), 2, f(:), 'filled')
colorbar
xlim([lon(eventID) - limitMesh, lon(eventID) + limitMesh])
ylim([lat(eventID) - limitMesh, lat(eventID) + limitMesh])
title(['Event ', num2str(eventID), ', Mw = ', num2str(Catalog.mag(eventID)), ...
        ', Model Run: ', cellStrRuns{iRun}])

idx1 = Catalog.flag > 0 & Catalog.t > Catalog.t(eventID) & Catalog.t < Catalog.t(eventID)+30;
idx2 = Catalog.t(Catalog.flag>0) > Catalog.t(eventID) & Catalog.t(Catalog.flag>0) < Catalog.t(eventID)+30;
scatter(lon(idx1), lat(idx1), 2+10*prob_trig_matrix(eventID,idx2), 'k', 'filled')

plot(lon(eventID), lat(eventID), 'r*')
hold on
if startsWith(Method.spatialDesign, 'aniso')
    plot([rupStart(1), rupEnd(1)], [rupStart(2), rupEnd(2)], 'r')
end
xlabel(['Sum of trigger prob: ', num2str(sum(prob_trig_matrix(eventID,idx2)))])

