function [ doubletPercentages, ...
           doubletSampleSize ] = compute_doubletPercentages( doubletDegree, ...
                                                             magnIntervals, ...
                                                             nMagnitudes, ...
                                                             idxEvents, ...
                                                             DoubletCriteria, ...
                                                             eventDate, ...
                                                             eventMag, ...
                                                             isSynthetic )
       
    %% Initialize output variables
    doubletSampleSize   = zeros(nMagnitudes, 1);
    doubletPercentages  = zeros(nMagnitudes, 3);
    
    %% Loop over all magnitude intervals
    for iRange = 1:nMagnitudes-1
        
        % For historic catalogs, check completeness of magnitude range
        if isSynthetic
            isComplete      = true(size(idxEvents));
        else
            yearComplete    = DoubletCriteria.yearComplete(iRange);
            isComplete      = year( eventDate(idxEvents) ) >= yearComplete;
        end
        
        % Extract doublet partner data for current magnitude interval
        iDegrees = doubletDegree( isComplete & eventMag(idxEvents) == magnIntervals(iRange), : );
        
        % Count number of doublets and multiplets
        doubletSampleSize(iRange)       = size(iDegrees, 1);
        doubletPercentages(iRange,:)    = [ sum( iDegrees==0 ), sum( iDegrees==1 ), sum( iDegrees(:,1)>1 )] ...
                                            / doubletSampleSize(iRange);
                                    
%         doubletPercentages(iRange,:)   = [sum(iDegrees(:,1)==0), sum(iDegrees(:,1)==1), sum(iDegrees(:,1)>1)] ...
%                                         / doubletSampleSize(iRange);
%         Doublet_perc04(iRange,:)    = [sum(iDegrees(:,2)==0), sum(iDegrees(:,2)==1), sum(iDegrees(:,2)>1)] ...
%                                         / doubletSampleSize(iRange);
%         Doublet_perc10(iRange,:)    = [sum(iDegrees(:,3)==0), sum(iDegrees(:,3)==1), sum(iDegrees(:,3)>1)] ...
%                                         / doubletSampleSize(iRange);
                                    
    end
       
end