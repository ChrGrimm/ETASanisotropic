function tempVolatility = evaluate_temporalClustering( Catalog, isTarget, startDate )

    iYears          = year( startDate + Catalog.t(isTarget) );
    iMonths         = month( startDate + Catalog.t(isTarget) );
    timeStamp       = iYears + (iMonths-1)/12;
    edgesTime       = (min(timeStamp)-1/24:1/12:max(timeStamp)+1/24);
    hist_temp       = histcounts( timeStamp, edgesTime );
    tempVolatility  = std(hist_temp)/mean(hist_temp);
    
end