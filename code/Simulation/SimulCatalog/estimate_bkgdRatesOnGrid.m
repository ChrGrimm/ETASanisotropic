function BackgrDistr = estimate_bkgdRatesOnGrid( TargetSettings_simul, ExperimentSettings )

    isBackgrDistrHomogen = ExperimentSettings.isBackgrDistrHomogen;
    
    %% Load estimation results
    pathResultsFolder = ['..\..\..\..\Estimation\', ExperimentSettings.strFolder_etasResults];
    clear TargetSettings
    load([pathResultsFolder, '\Results_estimation.mat'], 'Inputs', 'Catalog', 'ModelSummary', 'paramETAS')
    TargetSettings_estim = Inputs.TargetSettings;

    %% Precompute background grid
    [ gridX, gridY ]    = create_gridXY4polygon( TargetSettings_simul.polygonXY );
    BackgrDistr         = struct;
    BackgrDistr.gridX   = gridX';
    BackgrDistr.gridY   = gridY';
    uniqueGridX         = unique(gridX);
    uniqueGridY         = unique(gridY);
    BackgrDistr.deltaX  = abs( uniqueGridX(2)-uniqueGridX(1) );
    BackgrDistr.deltaY  = abs( uniqueGridY(2)-uniqueGridY(1) );
    
    %% Check whether estimation results can be used to estimate an inhomogeneous background distribution
    if ~isBackgrDistrHomogen
        
        % If simulation polygon is not inside estimation polygon, then set isBackgrDistrHomogen = true
        polygonTarget_estim = TargetSettings_estim.polygonTarget;
        polygonTarget_simul = TargetSettings_simul.polygonTarget;
        isSubPolygon = all( inpolygon(polygonTarget_simul(:,1), polygonTarget_simul(:,2), ...
                                      polygonTarget_estim(:,1), polygonTarget_estim(:,2)) );
        if ~isSubPolygon
            isBackgrDistrHomogen = true;
        end
    end
    
    %% Compute distribution of background events over spatial grid
    if isBackgrDistrHomogen
        nGridPoints        = length(gridX);
        BackgrDistr.prob   = (1/nGridPoints):(1/nGridPoints):1;
    
    else
        % Extract events from catalog within time window
        Catalog             = Catalog( Catalog.flag>=0, : );
        bandwidth           = precompute_backgrSpatDistr( Catalog, 'degree', 'aniso' );

        % Gaussian function of background rates (Jalilian, 2019) 
        dGauss              = @(r2, sig) exp(-r2 ./(2 * sig.^2)) ./ (2 * pi * sig.^2);
        % Calculate distances of events and grid points (nEvents x nGridPoints)
        dist2               = (Catalog.x - gridX).^2 + (Catalog.y - gridY).^2;
        % Calculate background rate per grid point
        bkgdRates           = sum( dGauss(dist2, bandwidth( Catalog.flag >= 0 )) .* Catalog.backgrProb, 1)';
        BackgrDistr.gridX   = gridX';
        BackgrDistr.gridY   = gridY';
        BackgrDistr.prob    = cumsum(bkgdRates)/sum(bkgdRates);
    
    end
    
    %% Compute mean number of background events in simulated time-space-magnitude window
    nDays                   = TargetSettings_simul.tWindow_days(2);
    polygonSize             = TargetSettings_simul.polyArea;
    magnScaling             = exp(paramETAS(10) * (TargetSettings_estim.Mc-TargetSettings_simul.Mc));
    
    % Scale number of background events
    BackgrDistr.nBackgrMean = ModelSummary.mu_backgrPerDayAndDegree * nDays * polygonSize * magnScaling;

end
%% Old code pieces
    %     % Scale background rates with area of grid cell
    %     if strcmp(Method.spaceUnit, 'degree')
    %         areaGridP   = dX * dY;
    %     else
    %         areaGridP   = (111 * dX * cosd(reshape(gridY, dim_xy))) *  (111 * dY);
    %     end
    %     
    %     bkgd = bkgd .* areaGridP / diff(SimulSettings.timeWindow);