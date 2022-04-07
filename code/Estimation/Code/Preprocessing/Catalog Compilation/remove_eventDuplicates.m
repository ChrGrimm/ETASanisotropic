function Catalog = remove_eventDuplicates( Catalog )
%% Delete events that occurred exactly at the same time (keep strongest candidate)

    [uniqueDates, ~, idxDates]  = unique(Catalog.date);
    idxSameTime                 = find(hist(idxDates, 1:length(uniqueDates))>1);
    idx2delete                  = [];
    for i=1:length(idxSameTime)
       idxEvents                = find(idxDates==idxSameTime(i));
       [~, idxMaxMagn]          = max(Catalog.mag(idxEvents));
       idxEvents(idxMaxMagn)    = [];
       idx2delete               = [idx2delete; idxEvents];
    end
    Catalog(idx2delete,:)       = [];
    
end