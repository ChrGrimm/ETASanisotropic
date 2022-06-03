function cdfPlot_numberOfAftershocks( strLegend, iRun, isLastModel )

    hold on
    
    %% Preprocessing
    % Load simulation results
    load('Simulation_results.mat', 'tAnalysis_synthCat', 'tAnalysis_observed')
    
    % Extract number of synthetic catalogs  
    nRealiz     = size(tAnalysis_synthCat, 1);
    
    % Extract number of events in the observed catalog
    nObserved   = tAnalysis_observed.nTargetEvents(1);
    
    % Define colors and line style for iterative plots (edit manually)
    [ color, lineStyle ] = define_colors4plotsInLoops( iRun );
    
    %% Display
    disp([strLegend, ': ', num2str(sum(tAnalysis_synthCat.nTargetEvents >= nObserved)/nRealiz*100), ...
            '\% of the synthetic catalogs exceed the observed number of events.'])
    
    %% Plotting
    % Plot cdf of number of simulated aftershocks (in target time-space window)
    pPred = plot( sort(tAnalysis_synthCat.nTargetEvents, 'ascend'), (1/nRealiz):(1/nRealiz):1, ...
                'Color', color, 'LineStyle', lineStyle, 'LineWidth', 1.5, 'DisplayName', strLegend );
    
    % Grid and icon display style
    pPred.Annotation.LegendInformation.IconDisplayStyle = 'off';
    grid on
    
    % Plot vertical line for number of events in observed catalog
    if isLastModel
        plot(nObserved*ones(1,2), [0,1], 'Color', [0.5 0.5 0.5], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', ['Observed: ', num2str(nObserved)])
    end
    
    %% Figure design
    xlim([0, 7000])
    make_titleLeftCorner('(a)')
%     title('Ridgecrest (10 days)')

    xlabel('Number of aftershocks')
    ylabel('Predicted cdf')
    xticks(sort([0,1000,2000,4000,5000,6000,7000, nObserved]))

    legend show
    legend('Location', 'Southeast')
    
end