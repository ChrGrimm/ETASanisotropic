function plot_rateVsStrike(sampleStrikes, sampleEpiPos, forwardRateContrib_sum, forwardRateContrib_logl, lineStyle)

%     close all
    
    for iEvent=1:size(sampleStrikes,1)
        
        figure(iEvent)
        
        for iEpiPos=1:length(sampleEpiPos)
            
            color = define_colors4plotsInLoops( iEpiPos );
            
            [ iStrikes, idxSort ]   = sort( sampleStrikes(iEvent,:) );
            iRates_sum              = squeeze( forwardRateContrib_sum(iEvent, idxSort, iEpiPos) );
            iRates_logl             = squeeze( forwardRateContrib_logl(iEvent, idxSort, iEpiPos) );
            
            %% Plot strikes vs contributed forward rate
%             subplot(1,2,1)
            title('Rate (sum) Vs. Strike')
            hold on
            plot(iStrikes, iRates_sum, 'Color', color, 'LineWidth', 1.5, 'LineStyle', lineStyle)

%             subplot(1,2,2)
%             title('Rate (log-likelihood) Vs. Strike')
%             hold on
%             plot(iStrikes, iRates_logl, 'LineWidth', 1.5)
            
        end
        
%         subplot(1,2,1)
        legend(cellstr(num2str(sampleEpiPos')), 'Location', 'Southwest')
        xlabel('Strike Sample (°)')
        ylabel('Sum of trigger rate contribution')
        
%         subplot(1,2,2)
%         legend(cellstr(num2str(sampleEpiPos')))
                
    end

end

%             % Extract maximizing and original strike
%             strike4maxRate  = sampleStrikes( idxSort(1) );
%             idxEvent2plot   = Catalog.id == id2plot(iEvent);
%             strikeOriginal  = Catalog.strike( idxEvent2plot );
%             % Find events in time window
%             isInTime        = Catalog.t - Catalog.t( idxEvent2plot ) > 0 & Catalog.t - Catalog.t( idxEvent2plot ) <= timeWindow;
%             timeDiff        = Catalog.t(isInTime) - Catalog.t(idxEvent2plot);
