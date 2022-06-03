function plot_typeOfAftershock( legendNames, iRun )

    %% Plot settings
    color       = define_colors4plotsInLoops( iRun );
    strLegend   = legendNames{iRun};
    
    %% Load data
    load('Simulation_results.mat','tAnalysis_synthCat')
%     load('Simulation_resultsMainshocksOnly.mat','tAnalysis_synthCat') % TEMPORARY
    [ centerMagnitudes, ...
       resultsBath, ...
       resultsDoublet ] = preprocess_typeOfAftershockResults( tAnalysis_synthCat );
    
    %% Plotting
%     % Bath
%     plot(centerMagnitudes, resultsBath(:,1), 'r--')
%     hold on
%     plot(centerMagnitudes, resultsBath(:,2), 'b--')
%     plot(centerMagnitudes, resultsBath(:,3), 'g--')
    % Doublets
    
    plot(centerMagnitudes, resultsDoublet(:,1), 'Color', color, 'LineStyle', '-.', ...
            'LineWidth', 1.5, 'DisplayName', [strLegend, ': independent'])
    hold on
%     plot(centerMagnitudes, resultsDoublet(:,2), 'Color', color, 'LineStyle', '--', ...
%             'LineWidth', 1.5, 'DisplayName', [strLegend, ': unlinked, same cluster'])
%     plot(centerMagnitudes, resultsDoublet(:,3), 'Color', color, 'LineStyle', '-', ...
%             'LineWidth', 1.5, 'DisplayName', [strLegend, ': trigger-related'])
    plot(centerMagnitudes, resultsDoublet(:,2)+resultsDoublet(:,3), 'Color', color, 'LineStyle', '-', ...
            'LineWidth', 1.5, 'DisplayName', [strLegend, ': same cluster'])
       
    nFontSize = 12;
    legend show
%     legend('Location', 'Southwest', 'FontSize', nFontSize)
    legend('Location', 'Northwest', 'NumColumns', 2)
%     title('Bath Law')
    ylabel('Share of events', 'FontSize', nFontSize)
    xlabel('Magnitude of triggering event (Mw)', 'FontSize', nFontSize)
%     a = get(gca,'XTickLabel');  
%     set(gca,'XTickLabel',a,'fontSize',nFontSize)
%     a = get(gca,'YTickLabel');  
%     set(gca,'YTickLabel',a,'fontSize',nFontSize)
%     ylim([0,1.05])
    ylim([0,1.15])
%     ylim([0,1])
    xlim([5.8,9])
    
     make_titleLeftCorner( '(d)' );
    
end