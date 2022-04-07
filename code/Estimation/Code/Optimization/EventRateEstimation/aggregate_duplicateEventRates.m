function eventRateData = aggregate_duplicateEventRates( eventRateData, idx, strMode )

    if strcmp( strMode, 'col' )
        eventRateData(:,idx(1)) = sum(eventRateData(:,idx'),2);
        eventRateData(:,idx(2)) = [];
        
    elseif strcmp( strMode, 'row' )
        eventRateData(idx(1),:) = sum(eventRateData(idx,:));
        eventRateData(idx(2),:) = [];
        
    else
        error('Unknown mode!')
        
    end

end