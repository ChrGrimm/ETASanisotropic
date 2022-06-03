function pdfPlot_magnOfStrongestAftersh( strLegend, iRun, isLastModel )

    hold on
    
    %% Preprocessing
    % Load simulation results
    load('Simulation_results.mat', 'tAnalysis_synthCat', 'tAnalysis_observed')
    
    % Extract number of synthetic catalogs  
    nRealiz     = size(tAnalysis_synthCat, 1);
    
    % Extract number of events in the observed catalog
    magObserved   = tAnalysis_observed.largestMagn;
    
    % Define colors and line style for iterative plots (edit manually)
    [ color, lineStyle ] = define_colors4plotsInLoops( iRun );
    
    %% Display
    disp([strLegend, ': ', num2str(sum(tAnalysis_synthCat.largestMagn >= magObserved-0.05)/nRealiz*100), ...
            '\% of the synthetic catalogs exceed the observed largest aftershock magnitude.'])
    
    %% Plots
    % Plot kernel density of pdf of simulated largest aftershock
    % magnitudes
    [f,x] = ksdensity( tAnalysis_synthCat.largestMagn, 2.5:0.1:7.5, 'Support', 'positive');
    pPred = plot(x, f, 'Color', color, 'LineStyle', lineStyle, 'LineWidth', 1.5, 'DisplayName', strLegend);
    
    % Grid and icon display style
    pPred.Annotation.LegendInformation.IconDisplayStyle = 'off';
    grid on
    
    % Alternative: plot cdf
%     plot( sort(tAnalysis_synthCat.largestMagn, 'ascend'), (1/nRealiz):(1/nRealiz):1, ...
%         'Color', color, 'LineStyle', lineStyle, 'LineWidth', 1.5, 'DisplayName', strLegend )

    % Plot vertical line for largest aftershock magnitude in observed catalog
    if isLastModel
        plot(magObserved*ones(1,2), [0,1], 'Color', [0.5 0.5 0.5], 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', ['Observed: ', num2str(magObserved)])
    end
    
    %% Figure design
    xlim([3, 7.5])
    make_titleLeftCorner('(b)')
%     title('Ridgecrest (10 days)')

    xlabel('Largest aftershock magnitude (Mw)')
    ylabel('Predicted pdf')
    xticks(unique(sort([2.5:0.5:7.5, magObserved])))
    
    legend show
    legend('Location', 'Northwest')
    
end