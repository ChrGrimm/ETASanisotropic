function [ doubletPercentages, ...
            nIntervalsQuantiles ] = group_doubletPerc4Catalogs( tAnalysis, ...
                                                                magnRanges2analyze, ...
                                                                magnIntervals, ...
                                                                intervals2Merge )

    %% Initialize auxiliary matrices and output matrix
    % Matrix dimensions
    nCatalogs           = size(tAnalysis,1);
    nIntervals          = length(magnIntervals)-1;
    nEvents             = zeros( nCatalogs, nIntervals );
    doubletPercentages  = zeros( nCatalogs, nIntervals );
    
    %% Extract original analysis data from input variables
    rawPercentages      = tAnalysis.percDoublets04;
    rawNumberEvents     = tAnalysis.nEventsPerMagnRange;
    rawIntervals        = magnRanges2analyze;
    
    %% Compute doublet percentages
    % Loop over realizations
    for iRealiz = 1:nCatalogs
        nDoublets = rawPercentages{iRealiz} .* rawNumberEvents{iRealiz};
        % Loop over magnitude intervals
        for jMagn = 1:nIntervals
            % Find analysis magnitudes that belong to current magnitude interval
            isInInterval = rawIntervals >= magnIntervals(jMagn) & ...
                            rawIntervals < magnIntervals(jMagn+1);
            % Compute aggregated doublet percentage in magnitude interval
            nEvents(iRealiz,jMagn)              = sum(rawNumberEvents{iRealiz}(isInInterval));
            doubletPercentages(iRealiz,jMagn)   = sum( sum( nDoublets(isInInterval,2:3), 'omitnan' ), 'omitnan' ) ...
                                                    / nEvents(iRealiz,jMagn);
        end
    end
    
    %% Combine magnitude ranges that belong together
    if ~isempty(intervals2Merge)
        idx                             = intervals2Merge;
        doubletPercentages(:,idx(1))    = sum( doubletPercentages(:,idx).*nEvents(:,idx), 2 ,'omitnan' ) ./ ...
                                            sum( nEvents(:,idx), 2 ,'omitnan' );
        doubletPercentages(:,idx(2))    = [];
        nIntervalsQuantiles             = size(doubletPercentages,2)-1;
    else
        nIntervalsQuantiles             = nIntervals;
    end
                                                         
end

