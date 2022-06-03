function [ doubletPerc_grouped, ...
            nIntervalsQuantiles ] = group_doubletPerc4Catalogs1( doubletPercentage, ...
                                                                 doubletSampleSize, ...
                                                                 doubletMagnIntervals, ...
                                                                 magnIntervals_grouped, ...
                                                                 idx_magnIntervals2merge )

    %% Initialize auxiliary matrices and output matrix
    % Matrix dimensions
    nCatalogs               = length(doubletPercentage);
    nIntervals              = length(magnIntervals_grouped)-1;
    nEvents                 = zeros( nCatalogs, nIntervals );
    doubletPerc_grouped     = zeros( nCatalogs, nIntervals );
    
    %% Compute grouped doublet percentages
    % Loop over realizations
    for iRealiz = 1:nCatalogs
        % Compute absolute number of doublets (percentage * sample size)
        doubletsNumber = doubletPercentage{iRealiz} .* doubletSampleSize{iRealiz};
        
        % Loop over magnitude intervals
        for jMagn = 1:nIntervals
            
            % Find analysis magnitudes that belong to current magnitude interval
            isInGroupedInterval = doubletMagnIntervals >= magnIntervals_grouped(jMagn) & ...
                                  doubletMagnIntervals <  magnIntervals_grouped(jMagn+1);
                       
            % Compute aggregated number of events in grouped magnitude interval
            nEvents(iRealiz,jMagn)              = sum( doubletSampleSize{iRealiz}(isInGroupedInterval) );
            
            % Compute aggregated doublet percentage in grouped magnitude interval
            doubletPerc_grouped(iRealiz,jMagn)  = sum( sum( doubletsNumber(isInGroupedInterval,2:3), 'omitnan' ), 'omitnan' ) ...
                                                    / nEvents(iRealiz,jMagn);
                                                
        end
        
    end
    
    %% Optional: Merge magnitude intervals
    if isempty(idx_magnIntervals2merge)
        
        % no merge
        nIntervalsQuantiles             = nIntervals;
        
    else
        
        % Merge intervals with following indices
        idx                             = idx_magnIntervals2merge;
        
        % Compute aggregated doublet percentage
        doubletPerc_grouped(:,idx(1))   = sum( doubletPerc_grouped(:,idx).*nEvents(:,idx), 2 ,'omitnan' ) ./ ...
                                            sum( nEvents(:,idx), 2 ,'omitnan' );
                                        
        % Remove additional column
        doubletPerc_grouped(:,idx(2))   = [];
        
        % Compute number of magnitude interval quantiles 
        nIntervalsQuantiles             = size(doubletPerc_grouped,2)-1;
        
    end  
                                                         
end

