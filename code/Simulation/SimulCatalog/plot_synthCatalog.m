function plot_synthCatalog( Catalog, ...
                            YearStat, ...
                            MagnStat, ...
                            SpaceStat, ...
                            edges_year, ...
                            edges_magn, ...
                            gridLon, ...
                            gridLat, ...
                            TargetWindow, ...
                            strBeta )    

    close all
    
    %% Prepare original catalog
    % Add year and absolute lon/lat information to synthetic catalog
    Catalog.evYears = year( TargetWindow.firstEvent + Catalog.t );

    % Extract events inside polygon boundaries
    isInPoly        = inpolygon( Catalog.x, Catalog.y, ...
                                 TargetWindow.polygonXY(:,1), TargetWindow.polygonXY(:,2) );
%     
    Catalog         = Catalog( isInPoly, : );
    isInTime        = Catalog.t >= TargetWindow.timeWindow(1) & Catalog.t <= TargetWindow.timeWindow(2);
    isHistTrig      = Catalog.backgrProb <= 0.5;
    
    %% Temporal occurrence statistic
    % Figure a): Common statistics for bkgd and trig
    % Calculate mean and 5/95% quantiles of YearStat
    yearStat_simu   = quantile(YearStat(:,:,1)+YearStat(:,:,2), [0.05, 0.5, 0.95], 2);
    edges_yearAll   = [min(Catalog.evYears):(edges_year(1)-1), edges_year];
    nEdges_pre      = length(edges_yearAll)-length(edges_year);
    
    yearStat_hist(:,1) = histcounts(Catalog.evYears(~isHistTrig), edges_yearAll)';    
    yearStat_hist(:,2) = histcounts(Catalog.evYears(isHistTrig), edges_yearAll)';     
       
    figure(1);
    bar(edges_yearAll(2:end)-1, yearStat_hist, 'stacked', 'LineStyle', 'none')
    hold on
    plot(edges_year(2:end)-1, yearStat_simu(:,2), 'k-', 'LineWidth', 1.5)
    plot(edges_year(2:end)-1, yearStat_simu(:,1), 'k--', 'LineWidth', 1.5)
    plot(edges_year(2:end)-1, yearStat_simu(:,3), 'k--', 'LineWidth', 1.5)
    legend('historic: P(background)>0.5', 'historic: P(triggered)>=0.5', ...
            'synthetic: all (median)', 'synthetic: all (5%/95%)')
    title({'Number of event occurrences in 2-year windows','Synthetic Catalog from 1966/1/1'})
    % title('Number of event occurrences in 2-year windows: Synthetic Catalog from 1966/1/1')
    xlabel('Time (years)')
    ylabel('Number of occurrences')
    axis([min(Catalog.evYears)-1, max(Catalog.evYears)+1, 0, 800])

    saveas(gcf, ['nOcc_', strBeta, '.png']);
    
    % Figure b): Separate statistics for bkgd and trig
    yearStat_quanBkgd   = quantile(YearStat(:,:,1), [0.05, 0.5, 0.95], 2);
    yearStat_quanTrig   = quantile(YearStat(:,:,2), [0.05, 0.5, 0.95], 2);
    
    figure(2);
    subplot(1,2,1)
    plot(edges_year(2:end)-1, yearStat_hist((1+nEdges_pre):end,1), 'bx-', 'LineWidth', 1.5)
    hold on
    plot(edges_year(2:end)-1, yearStat_quanBkgd(:,2), 'r-', 'LineWidth', 1.5)
    plot(edges_year(2:end)-1, yearStat_quanBkgd(:,1), 'r--', 'LineWidth', 1.5)
    plot(edges_year(2:end)-1, yearStat_quanBkgd(:,3), 'r--', 'LineWidth', 1.5)
    legend('historic: P(background)>0.5', 'synthetic: background (median)', 'synthetic: background (5%/95%)')
    title({'Number of background event occurrences in 2-year windows','Synthetic Catalog from 1966/1/1'})
    xlabel('Time (years)')
    ylabel('Number of occurrences')
    axis([edges_year(1)-2, edges_year(end)+2, 0, 500])

    saveas(gcf, ['nTrigBkgd_', strBeta, '.png']);
    
    subplot(1,2,2)
    plot(edges_year(2:end)-1, yearStat_hist((1+nEdges_pre):end,2), 'bx-', 'LineWidth', 1.5)
    hold on
    plot(edges_year(2:end)-1, yearStat_quanTrig(:,2), 'r-', 'LineWidth', 1.5)
    plot(edges_year(2:end)-1, yearStat_quanTrig(:,1), 'r--', 'LineWidth', 1.5)       
    plot(edges_year(2:end)-1, yearStat_quanTrig(:,3), 'r--', 'LineWidth', 1.5)
    legend('historic: P(triggered)>0.5', 'synthetic: triggered (median)', 'synthetic: triggered (5%/95%)')
    title({'Number of triggered event occurrences in 2-year windows','Synthetic Catalog from 1966/1/1'})
    xlabel('Time (years)')
    ylabel('Number of occurrences')
    axis([edges_year(1)-2, edges_year(end)+2, 0, 500])
    
    %% Plot single simulation outputs
%     figure(3);
%     
%     for i=1:8
%         subplot(3,3,i);
%         data = [ YearStat(:,i,1), YearStat(:,i,2) ];
%         bar(edges_yearAll(2:nEdges_pre+1)-1, sum(yearStat_hist(1:nEdges_pre,:), 2), 'stacked', 'LineStyle', 'none')
%         hold on
%         bar(edges_year(2:end)-1, data, 'stacked', 'LineStyle', 'none') 
%         plot(edges_year(2:end)-1, sum(yearStat_hist(1+nEdges_pre:end,:), 2), 'k-', 'LineWidth', 1.5)
%         xticks([1965, 1985, 2005])
%                
%     end
%     subplot(3,3,9);
% %     title({'Number of event occurrences in 2-year windows','Exemplary synthetic catalogs from 1966-2007'})
%     data = [ YearStat(:,i,1), YearStat(:,i,2) ];
%     bar(edges_yearAll(2:nEdges_pre+1)-1, sum(yearStat_hist(1:nEdges_pre,:), 2), 'stacked', 'LineStyle', 'none')
%     hold on
%     bar(edges_year(2:end)-1, data, 'stacked', 'LineStyle', 'none') 
%     plot(edges_year(2:end)-1, sum(yearStat_hist(1+nEdges_pre:end,:), 2), 'k-', 'LineWidth', 1.5)
%     axis([0, 1, 0, 1])
%     axis off
%     legend('historic: 1960-1965', 'synthetic: backgr.', 'synthetic: trigg.', 'historic: 1966-2007', 'Location', 'Northwest')

    figure(3);
    
    for i=1:3
        subplot(2,2,i);
        data = [ YearStat(:,i,1), YearStat(:,i,2) ];
        bar(edges_yearAll(2:nEdges_pre+1)-1, [0 0 0], 'stacked', 'LineStyle', 'none')
        hold on
        bar(edges_year(2:end)-1, data, 'stacked', 'LineStyle', 'none') 
       
        plot(edges_year(2:end)-1, sum(yearStat_hist(1+nEdges_pre:end,:), 2), 'k-', 'LineWidth', 1.5)
        xticks([1965, 1985, 2005])
        xlabel('Year')
        ylabel('Number of occurrences')
               
    end
    subplot(2,2,4);
%     title({'Number of event occurrences in 2-year windows','Exemplary synthetic catalogs from 1966-2007'})
    data = [ YearStat(:,i,1), YearStat(:,i,2) ];
    bar(edges_yearAll(2:nEdges_pre+1)-1, sum(yearStat_hist(1:nEdges_pre,:), 2), 'stacked', 'LineStyle', 'none')
    hold on 
    L1 = bar(edges_year(2:end)-1, data, 'stacked', 'LineStyle', 'none');     
    L2 = plot(edges_year(2:end)-1, sum(yearStat_hist(1+nEdges_pre:end,:), 2), 'k-', 'LineWidth', 1.5);
    axis([0, 1, 0, 1])
    axis off
    legend([L1, L2], 'synthetic: backgr.', 'synthetic: trigg.', 'historic: 1973-2019', 'Location', 'Northwest')


    saveas(gcf, ['setSynthCats_', strBeta, '.png']);

    %% Magnitude statistic
    % Cumulate synthetic catalog statistics
    MagnStat_exc        = flipud(cumsum(MagnStat(end:-1:1,:), 1));
    % Cumulate historic catalog statistics
    MagnStat_hist       = histcounts(Catalog.mag(isInTime)+4.5, edges_magn)';
%     MagnStat_bkgd       = histcounts(Catalog.mag(~isHistTrig&isInTime)+4.5, edges_magn)';
    MagnStat_excHist    = flipud(cumsum(MagnStat_hist(end:-1:1)));
%     MagnStat_excBkgd    = flipud(cumsum(MagnStat_bkgd(end:-1:1)));

    % Calculate mean and 5/95% quantiles of YearStat
    MagnStat_quant      = quantile(MagnStat_exc(:,:), [0.05, 0.5, 0.95], 2);
    % Compute log10 of data
    MagnStat_quant      = log10(MagnStat_quant);
    MagnStat_excHist    = log10(MagnStat_excHist);
%     MagnStat_excBkgd    = log10(MagnStat_excBkgd);
%     % Scale 
%     MagnStat_quant      = MagnStat_quant./max(MagnStat_quant);
%     MagnStat_excTrig    = MagnStat_excTrig/max(MagnStat_excTrig);
%     MagnStat_excBkgd    = MagnStat_excBkgd/max(MagnStat_excBkgd);

    
    MagnStat_quant(isinf(MagnStat_quant))       = NaN;
    MagnStat_excHist(isinf(MagnStat_excHist))   = NaN;
%     MagnStat_excBkgd(isinf(MagnStat_excBkgd))   = NaN;
    
    figure(4);
    plot(edges_magn(2:end)-0.05, MagnStat_excHist, 'blue', 'LineWidth', 1.5)
    hold on
%     plot(edges_magn(2:end)-0.05, MagnStat_excTrig, 'red', 'LineWidth', 1.5)
    plot(edges_magn(2:end)-0.05, MagnStat_quant(:,2), 'k-', 'LineWidth', 1.5)
    plot(edges_magn(2:end)-0.05, MagnStat_quant(:,1), 'k--', 'LineWidth', 1.5)
    plot(edges_magn(2:end)-0.05, MagnStat_quant(:,3), 'k--', 'LineWidth', 1.5)
    legend('historical', 'synthetic: median', 'synthetic: 5/95%')
    title({'Log10 of Magnitude Frequency Distribution (Exceedance)','Historical: 1926-2007; Synthetic: 1966-2007'})
    % title('Number of event occurrences in 2-year windows: Synthetic Catalog from 1966/1/1')
    xlabel('Magnitudes (Mw)')
    ylabel('Log10(Number of occurrences)')

    saveas(gcf, ['magnDistr_', strBeta, '.png']);
 
    
    %% Space statistics    
    SpaceStat_bkgd = mean(SpaceStat(:,:,1),2);
    SpaceStat_trig = mean(SpaceStat(:,:,2),2);
    gridX               = gridLon(:)';
    gridY               = gridLat(:)'; 
    % bkgd
    dist2               = (Catalog.x - gridX).^2 + (Catalog.y - gridY).^2;
    [~, gridIdx]        = min(dist2,[], 2);
    OrigSpaceStat       = histcounts(gridIdx, 0.5:1:length(gridX)+0.5)';

    % Back-convert event locations to regular lon-lat  
    cntLon = 136.5;
    cntLat = 36;
    gridLat1 = gridLat + cntLat;
    gridLon1 = gridLon / cos(cntLat * pi / 180) + cntLon;
    % Background
    figure;
    scatter(gridLon1(:), gridLat1(:), 5, SpaceStat_bkgd, 'filled');
    hold on
    plot_shapeMap
    axis([128, 145, 27, 45])
    title('Mean number of events per grid cell: Synthetic background events')
    colorbar
    % Triggered
    figure;
    scatter(gridLon1(:), gridLat1(:), 5, SpaceStat_trig, 'filled');
    hold on
    plot_shapeMap
    axis([128, 145, 27, 45])
    title('Mean number of events per grid cell: Synthetic triggered events')
    colorbar

    % Original catalog
    figure;
    scatter(gridLon1(:), gridLat1(:), 5, log10(OrigSpaceStat), 'filled');
    hold on
    plot_shapeMap
    axis([128, 145, 27, 45])
    xlabel('Longitude')
    ylabel('Latitude')
    title('Log10 of event occurrences per grid cell: Original catalog')
    colorbar
    lim = caxis;

    % Both
    figure;
    i = round(100*rand);
    SpaceStat_i = SpaceStat(:,i,1) + SpaceStat(:,i,2);
%     scatter(gridLon1(:), gridLat1(:), 5, SpaceStat_trig+SpaceStat_bkgd, 'filled');
    scatter(gridLon1(:), gridLat1(:), 5, log10(SpaceStat_i), 'filled');
    hold on
    plot_shapeMap
    axis([128, 145, 27, 45])
    xlabel('Longitude')
    ylabel('Latitude')
    title('Log10 of event occurrences per grid cell: Exempl. synthetic event set')
    colorbar
    caxis(lim)
