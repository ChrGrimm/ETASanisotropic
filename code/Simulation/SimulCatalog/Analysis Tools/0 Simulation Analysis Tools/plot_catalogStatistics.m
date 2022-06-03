function plot_catalogStatistics( SynthCatStatistics, strModelNames )

    %% Settings
    strVersion          = 'poster';
    strData             = 'tempVolatility';
%     strData             = 'nTargetEvents';
%     strData             = 'largestMagn';    
%     strData             = 'spatVolatility';
    
    %% Process synthetic catalog statistics
    % Extract desired column from synthetic catalog statistics
    dataSimulation      = SynthCatStatistics.(strData);
    % Extract number of models
    nModels             = size(SynthCatStatistics.nTargetEvents, 2);
    
    %% Load historic catalog results
    load('Results_histCatalogs.mat', 'tAnalysis_historicCat')
    % Extract results for desired historic catalog
    idHistCat           = 1;
    dataHistoric        = tAnalysis_historicCat.(strData);
    dataHistoric        = dataHistoric(idHistCat);
    
    %% Plotting
    figure
    % Boxplot of simulation data
    boxplot(dataSimulation, 1:nModels) 
    hold on
%     plot(1:nModels, CatStatistics.integrRateTempVol,'b*','MarkerSize',14,'DisplayName','Integrated Rate')

    %% Line plot of historic catalog
    plot([0.5:4.5], dataHistoric(1) * ones(1,5),'k-', 'LineWidth', 1.5, 'DisplayName','orig. catalog')
%     plot([4.5:8.5], dataHistoric(2) * ones(1,5),'k--','DisplayName','orig. catalog (CAL)')
    
    %% Plot settings
    if strcmp(strVersion, 'paper')
        
%         legend show
%         legend('Location', 'Northwest')
        ylim=get(gca,'ylim');
%         text(0.9, ylim(2)-0.15, 'CAL', 'FontSize', 12) %, 'FontWeight', 'bold')
        % text(0.9, ylim(2)-0.15, 'Southern California', 'FontSize', 10) %, 'FontWeight', 'bold')
%         xticklabels({'M_0', 'M_1', 'M_2', 'M_3'})
        xticklabels({'Poisson', 'Standard ETAS', 'Modified ETAS'})
        yticks(0:0.5:ylim(2))
%       title('Box Plot: Temporal volatility (coefficient of variance)')
        t = title('(a)', 'FontSize', 10);
        t.Units = 'Normalize'; 
        t.Position(1) = 0;
        t.HorizontalAlignment = 'left';
%     ylim([-inf, 1.1*tAnalysis_historicCat.tempVolatility(idOrigCat)])

    else
        
        nFontsize = 24;
        legend('Location', 'Southoutside', 'Orientation', 'horizontal', 'FontSize', nFontsize)
%         legend show
%         legend('Location', 'Southoutside', 'FontSize', nFontsize)
        title('Temporal clustering', 'FontSize', nFontsize)
        yticks([0.5:0.5:2])
        xticklabels({'M_0', 'M_1', 'M_2', 'M_3'})
%         xticklabels({'PSHA', 'ETAS_0', 'ETAS_1'})
        ylabel('Coefficient of Variation', 'FontSize', nFontsize)
        a = get(gca,'XTickLabel');  
        set(gca,'XTickLabel',a,'fontsize',nFontsize)
        a = get(gca,'YTickLabel');  
        set(gca,'YTickLabel',a,'fontsize',nFontsize)
        
    end

end