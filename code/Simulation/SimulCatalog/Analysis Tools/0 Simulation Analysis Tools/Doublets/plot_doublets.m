function plot_doublets( legendNames, iRun )

    %% Specific settings
    % Plot version
    strVersion  = 'single plot'; % 'shared plot';
    % Investigated region
    load('Estimation_results.mat','TargetWindow')
    % Extract legend entry
    strLegend   = legendNames{ iRun };
    
    %% Define plot settings
    [ magnIntervals_grouped, ...
      idx_magnIntervals2merge, ...
      strXtickLabels, ...
      idHistCat ] = define_plotSettings4doublets( TargetWindow.region, strVersion );
    % Line color
    color       = define_colors4plotsInLoops( iRun );
    
    %% Load data
    % Simulation results
    load('Simulation_results.mat','tAnalysis_synthCat','tAnalysis_synthSeq')
%     load('Simulation_resultsMainshocksOnly.mat','tAnalysis_synthCat','tAnalysis_synthSeq') % TEMPORARY
%     tAnalysis_synthSeq = tAnalysis_synthSeq365;
    % Historic results
    addpath(genpath('C:\Work\PhD_Christian\Matlab\Results\Historic Catalogs'))
    load('Results_histCats_mainshockOnly.mat', 'tAnalysis_historicCat', 'HistoricCatalogs') % TEMPORARY
%     load('Results_histCatalogs.mat', 'tAnalysis_historicCat', 'HistoricCatalogs', 'DoubletCriteria')
    % Doublet criteria
    load('Simulation_settings.mat','DoubletCriteria')
%     DoubletCriteria.magnIntervals = DoubletCriteria.magnRanges2analyze; % TEMPORARY
    
    %% Preprocess data
    % Extract relevant historic catalog results                                                    
    tAnalysis_historicCat = tAnalysis_historicCat( idHistCat, : );
    
    % Evaluate doublet percentages in chosen magnitude intervals
    [ doubletPerc_simulation, ...
      doubletPerc_sequence, ...
      doubletPerc_historic, ...
      magnIntervals, ...
      nIntervalsQuantiles ] = preprocess_doubletPercentages4plot( tAnalysis_synthCat, ...
                                                                  tAnalysis_synthSeq, ...
                                                                  tAnalysis_historicCat, ...
                                                                  DoubletCriteria.magnIntervals(1:end-1), ...
                                                                  magnIntervals_grouped, ...
                                                                  idx_magnIntervals2merge );
    
    %% Plotting 
    if strcmp(strVersion, 'single plot')
        
        figure
        % Simulation results
%         subplot(1,2,iRun)
        hold on
        plot( [1:nIntervalsQuantiles, nIntervalsQuantiles+0.5], mean(doubletPerc_simulation, 'omitnan'), ...
                'color', 'b', 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', [strLegend, ' simulations'] )
        scatter( [1:nIntervalsQuantiles, nIntervalsQuantiles+0.5], mean(doubletPerc_simulation, 'omitnan'), ...
                25, 'b', 'filled', 'HandleVisibility', 'off' )
        % Plot confidence intervals
        fill( [1:nIntervalsQuantiles, nIntervalsQuantiles:-1:1], ...
              [zeros(1,nIntervalsQuantiles), fliplr(quantile(doubletPerc_simulation(:,1:nIntervalsQuantiles), 0.90))], ...
              'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'DisplayName', [strLegend, ' simulations 10/90%'])
          
        % Historic results            
        LineStyle = {'-', '--', ':'};
        displayNames = {'Original catalog', 'ISC-GEM (Japan)', 'ISC-GEM (Global)'};
        for iOrig = 1:length(idHistCat)
            plot( 1:size(doubletPerc_historic,2), doubletPerc_historic(iOrig,:), ...
            'color', 'k', 'LineStyle', LineStyle{iOrig}, ...
            'DisplayName', displayNames{iOrig} )
            scatter( 1:size(doubletPerc_historic,2), doubletPerc_historic(iOrig,:), ...
                     15, 'k', 'filled', 'HandleVisibility', 'off' )
            %'color', 'k', 'LineStyle', 'none', 'Marker', 'o', 'HandleVisibility', 'off' )
        end
        
        ylim([0,0.4])
%         ylim([0,0.75]) 

        strVersion = 'paper';
        if strcmp(strVersion, 'paper')
            % Legend and title
            legend show
            legend('Location', 'Northeast')
            ylabel('Doublet percentage')
            xlabel('Magnitude of triggering event (Mw)')
            xticks([1:nIntervalsQuantiles, nIntervalsQuantiles+0.5])
            xticklabels(strXtickLabels)
            xlim([1,size(doubletPerc_simulation,2)])
            
            make_titleLeftCorner( '(c)' )
            
        elseif strcmp(strVersion, 'poster')          
            nFontSize = 24;
            title('M_0: Doublet Frequencies', 'FontSize', nFontSize)
            ylim([0,0.45]) % ylim([0,0.65])
            xlim([1,size(doubletPerc_simulation,2)-1])  
            xticks([1:nIntervalsQuantiles])
            xticklabels(strXtickLabels(1:end-1))
            ylabel('Doublet percentage', 'FontSize', nFontSize)
%             xlabel('Magnitude of triggering event (Mw)', 'FontSize', nFontSize)  
            xlabel('Magnitude (Mw)', 'FontSize', nFontSize)     
            a = get(gca,'XTickLabel');  
            set(gca,'XTickLabel',a,'fontsize',nFontSize)
            a = get(gca,'YTickLabel');  
            set(gca,'YTickLabel',a,'fontsize',nFontSize)
%             legend show
%             legend('Location', 'Northeast')
%             if iRun==2
% %                 Legend and title
%                 legend show
%                 legend('Location', 'Southoutside', 'Orientation', 'horizontal', 'FontSize', nFontSize)                
%             end
        end
        
    elseif strcmp(strVersion, 'shared plot')
        
        % Synthetic sequence results
        plot( (magnIntervals(1:end-1)+magnIntervals(2:end))/2, doubletPerc_sequence, ...
            'color', color, 'LineStyle', '--', 'LineWidth', 1.5, 'DisplayName', [strLegend, ' sequences'] )
        hold on
        % Synthetic catalog results
        plot( (magnIntervals(1:end-1)+magnIntervals(2:end))/2, mean(doubletPerc_simulation, 'omitnan'), ...
                'color', color, 'LineWidth', 1.5, 'DisplayName', [strLegend, ' catalogs'] )
        scatter( (magnIntervals(1:end-1)+magnIntervals(2:end))/2, mean(doubletPerc_simulation, 'omitnan'), ...
                15, color, 'filled', 'HandleVisibility', 'off' )
        % Legend and title
        legend show
        legend('Location', 'Northeast', 'FontSize', 9)
    %     title('Bath Law')
        ylabel('Doublet percentage')
        xlabel('Magnitude of triggering event (Mw)')
        if strcmp(TargetWindow.region, 'JPN')
            xlim([5.9,9.0])
        else strcmp(TargetWindow.region, 'CAL')
            xlim([5.9,7.5])
        end
        ylim([0, 0.18])
        
        make_titleLeftCorner( '(a)' )
        
    end
    
end