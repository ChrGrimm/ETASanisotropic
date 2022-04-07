% Catalog1 = Catalog(ismember(Catalog.id, TriggerRelationsETAS.evID(Catalog.mag(Catalog.flag>0)>=5.5 & TriggerRelationsETAS.nOffsprings_sumProb==0)),:);
% Catalog2 = Catalog1(Catalog1.depth < 80, :);

Catalog2 = Catalog(Catalog.mag>=6.5,:);

for i=1:size(Catalog2,1)
    idx                         = find(Catalog.id==Catalog2.id(i));
    rupLength                   = estimate_rupSize(Catalog.mag(idx), 'subduction', 'all', 'degree');
    idxInTime                   = find(Catalog.date > Catalog.date(idx) & Catalog.date < Catalog.date(idx)+days(1));
    distances                   = getDistance(Catalog.lon(idxInTime), Catalog.lat(idxInTime), ...
                                              Catalog.lon(idx), Catalog.lat(idx), 'degree');
    xTimesRupL                  = distances / rupLength;
    isInSpace                   = xTimesRupL < 10;
    idxInTimeSpace              = idxInTime(isInSpace);
    distances                   = distances(isInSpace);
    
    % Create output table
    NearestEvents               = table;
    [ NearestEvents.dist, ...
        idxSort ]               = sort(distances);
    NearestEvents.xRestrFactor  = NearestEvents.dist / rupLength;
    tDifferences                = days(Catalog.date(idxInTimeSpace) - Catalog.date(idx));
    NearestEvents.tDiff         = tDifferences(idxSort);

    TriggerEvent                = Catalog(idx,:);
    TriggerEvent.mag
    
end

% %% Plot
% plot_shapeMapscatter(Catalog1.lon, Catalog1.lat, 20, 'b', 'filled')
% scatter(Catalog2.lon, Catalog2.lat, 20, 'r', 'filled')
% scatter(Catalog.lon(Catalog.mag>=5.5), Catalog.lat(Catalog.mag>=5.5), 5, 'k', 'filled')
% scatter(Catalog.lon(Catalog.mag>=6.5), Catalog.lat(Catalog.mag>=6.5), 15, 'g', 'filled')
% load('Estimation_results.mat', 'Inputs')
% polygonTarget = Inputs.TargetSettings.polygonTarget;
% polygonComple = Inputs.TargetSettings.polygonComple;
% plot(polygonTarget(:,1), polygonTarget(:,2), 'k-')
% 
% Catalog55 = Catalog(Catalog.mag>=5.5, :);
% 
% Trigger55 = table;
% Trigger55.id = Catalog55.id(Catalog55.flag>0);
% Trigger55.mag = Catalog55.mag(Catalog55.flag>0);
% Trigger55.depth = Catalog55.depth(Catalog55.flag>0);
% Trigger55.nOffsprings = TriggerRelationsETAS.nOffsprings_sumProb(ismember(TriggerRelationsETAS.evID, Catalog55.id(Catalog55.flag>0)));
% figure
% scatter(Trigger55.depth, Trigger55.nOffsprings, 10, Trigger55.mag, 'filled'); colorbar