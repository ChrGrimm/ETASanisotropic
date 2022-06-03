function plot_doubletsOrigCatalogs( magnIntervals, intervBelongTogether, idOrigCatalogs )
% CURRENTLY NOT USED

    load('Results_histCatalogs.mat', 'tAnalysis_historicCat', 'HistoricCatalogs', 'DoubletCriteria')
    LineStyle = {'-', '--', ':'};
    
    if ~isempty(intervBelongTogether)
        magnIntervals = magnIntervals([1:3,intervBelongTogether(1),intervBelongTogether(2)+1]);
    end
    
    % Extract original analysis data from input variables
    rawPercentages      = tAnalysis_historicCat.percDoublets04;
    rawNumberEvents     = tAnalysis_historicCat.nEventsPerMagnRange;
    rawIntervals        = DoubletCriteria.magnRanges2analyze(1:end-1);
    
    for iOrig = 1:length(idOrigCat)
        
        plot( (1:4)', sum(percData{idOrigCatalogs(iOrig)}(:,2:3), 2), ...
              'Color', 'k', 'LineStyle', LineStyle{iOrig}, ...
              'DisplayName', HistoricCatalogs.catName{idOrigCat(iOrig)} )
        
    end
    
end