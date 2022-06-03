function scatterPlot_aftershLocations

    %% Preprocessing
    % Load simulation inputs and results
    load('Inputs_simulation.mat', 'SimulSettings')
    load('Results_simulation.mat', 'tAnalysis_synthCat', 'Observations')
    load('C:\Work\PhD_Christian\Projects\ETAS Model\Results\2nd Paper\211131 Submission\Loglik Space Fit\loglik0_Rid64_spaceM71_10days.mat', 'loglik0')
    load(SimulSettings.strPolygon, 'xGrid', 'yGrid')
%     load('C:\Work\PhD_Christian\Projects\ETAS Model\Results\2nd Paper\211131 Submission\ObsRidgecrestM71_10days.mat', 'Observations')
    
    % Extract number of synthetic catalogs and spatial grid cells
    nRealiz     = size(tAnalysis_synthCat, 1);
    nGridCells  = length(tAnalysis_synthCat.spatDistr{1});
    
    %% Compute averaged spatial aftershock distributions from n realizations
    % Compile matrix from n realizations
    SpatialDistr = zeros(nRealiz, nGridCells);
    for iRealiz = 1:nRealiz
        SpatialDistr(iRealiz,:) = tAnalysis_synthCat.spatDistr{iRealiz};
    end
    
    % Compute mean distribution
    meanSpatialDistr = mean(SpatialDistr)/sum(mean(SpatialDistr));
    
    %% Compute loglikelihood of spatial aftershock forecast
    nCounts     = round( tAnalysis_synthCat.spatDistr{1}*tAnalysis_synthCat.nTargetEvents(1) );
    isZero      = nCounts==0; %meanSpatialDistr==0;
    if find(meanSpatialDistr==0 & nCounts>0)
        warning('Observed event occurs at grid point that has P_j=0.')
        isZero = meanSpatialDistr==0 | nCounts==0;
    end

    loglik      = sum(nCounts(~isZero) .* log(meanSpatialDistr(~isZero)) ); % likelihood = prod(meanSpatialDistr.^nCounts);
    infoGain    = (loglik - loglik0) / tAnalysis_synthCat.nTargetEvents(1);
    disp(['Information Gain: ', num2str(infoGain)])
    
    %% Convert locations from centered x-y back to lat-lon scale
    centerLat   = Observations.lat(end);
    centerLon   = Observations.lon(end);
    latFactor  	= cos( centerLat * pi / 180 );
    
    % Convert grid locations
    lonGrid     = xGrid/latFactor+centerLon;
    latGrid     = yGrid + centerLat;
    
    % Convert observed event locations
    isFlag      = Observations.flag==1;
    lonObserved = Observations.x(isFlag)/latFactor + centerLon;
    latObserved = Observations.y(isFlag) + centerLat;
    
    %% Plotting
    figure
    a = axes;
    
    % Plot averaged spatial aftershock distribution on grid
    scatter(lonGrid, latGrid, 10, log10(meanSpatialDistr), 'filled')
    
    % Set limits of colorbar
    set(a, 'CLim', [-7.5 -1.5]);
    hold on
    colorbar

    % Plot observed aftershock locations
    scatter(lonObserved, latObserved, 1, 'k', 'filled')
    
    %% Figure design
    x_minmax = [min(lonGrid), max(lonGrid)] + [-0.025, 0.025];
    y_minmax = [min(latGrid), max(latGrid)] + [-0.025, 0.025];
    xlim(x_minmax)
    ylim(y_minmax)
    
    xticks(floor(xmin):0.2:ceil(xmax))
    yticks(floor(ymin):0.2:ceil(ymax))
    xlabel('Longitude (°)')
    ylabel('Latitude (°)')
    
    make_titleLeftCorner('(c)')
        
end