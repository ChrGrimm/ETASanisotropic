function [ doubletPerc_simulation, ...
           doubletPerc_sequence, ...
           doubletPerc_historic, ...
           magnIntervals_grouped, ...
           nIntervalsQuantiles ] = preprocess_doubletPercentages4plot( tAnalysis_synthCat, ...
                                                                       tAnalysis_synthSeq, ...
                                                                       tAnalysis_historicCat, ...
                                                                       doubletMagnIntervals, ...
                                                                       magnIntervals_grouped, ...
                                                                       idx_magnIntervals2merge )
    %% Evaluate doublet percentages in synthetic catalog results
%     [ doubletPerc_simulation, ...
%         nIntervalsQuantiles ] = group_doubletPerc4Catalogs( tAnalysis_synthCat, ...
%                                                             doubletMagnIntervals, ...
%                                                             magnIntervals_grouped, ...
%                                                             idx_magnIntervals2merge );
                                                        
    [ doubletPerc_simulation, ...
        nIntervalsQuantiles ] = group_doubletPerc4Catalogs1( tAnalysis_synthCat.doubletPercentage, ...
                                                             tAnalysis_synthCat.doubletSampleSize, ...
                                                             doubletMagnIntervals, ...
                                                             magnIntervals_grouped, ...
                                                             idx_magnIntervals2merge );
    
    %% Evaluate doublet percentages in synthetic sequence results
    doubletPerc_sequence = group_doubletPerc4Sequences( tAnalysis_synthSeq, magnIntervals_grouped );
                                                        
    %% Modify magnitude intervals if parts are to be merged
    if ~isempty(idx_magnIntervals2merge)
        magnIntervals_grouped = magnIntervals_grouped([1:3,idx_magnIntervals2merge(1),idx_magnIntervals2merge(2)+1]);
    end
    
    %% Evaluate doublet percentages for historic catalogs
    doubletPerc_historic = group_doubletPerc4Catalogs1( tAnalysis_historicCat.doubletPercentage, ...
                                                        tAnalysis_historicCat.doubletSampleSize, ...
                                                        doubletMagnIntervals, ...
                                                        magnIntervals_grouped, ...
                                                        [] );   
       
end