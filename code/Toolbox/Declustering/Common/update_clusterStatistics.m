function Clusters = update_clusterStatistics( Clusters, newCluster, Catalog, method )

    warning off
    
    %% Preprocessing
    % Row index of new cluster
    iC              = size(Clusters,1) + 1;
    
    % Cluster event indices
    idxEvents       = find( ismember(Catalog.id, newCluster.evIDs) );
    [~, idxSorted]  = sort(Catalog.mag(idxEvents),'descend');
    idxEventsSorted = idxEvents( idxSorted );
    
    % Join clusters (only for window method)
    if strcmp(method, 'window method') && ~isempty(Clusters)
        [iC, idxEventsSorted] = merge_clusters( Catalog, Clusters, iC, idxEventsSorted );
    end
    
    %% Update Cluster iC
    % General cluster information
    Clusters.clusterID(iC)          = iC;
    Clusters.eventIDs{iC}           = Catalog.id(idxEventsSorted);
    Clusters.size(iC)               = length(idxEventsSorted);

    % Mainshock information
    idxMainshock                    = idxEventsSorted(1);
    Clusters.mId(iC)           = Catalog.id(idxMainshock);
    Clusters.mDate(iC)         = Catalog.date(idxMainshock);
    Clusters.mMag(iC)         = Catalog.mag(idxMainshock);
%     Clusters.mLon(iC)          = Catalog.lon(idxMainshock);
%     Clusters.mLat(iC)          = Catalog.lat(idxMainshock);
%     Clusters.tWindow{iC}   = [ days(min(Catalog.date(idxEventsSorted))-Clusters.mainshDate(iC)), ...
%                                         days(max(Catalog.date(idxEventsSorted))-Clusters.mainshDate(iC)) ];
                                    
    % Strongest aftershock information
    idxAftershock                   = idxEventsSorted(2);
%     Clusters.aftershID(iC)      = Catalog.id(idxAftershock);
    Clusters.aMag(iC) = Catalog.mag(idxAftershock);
    
    %% Method-specific evaluations
    if strcmp(method, 'window method')
        Clusters.deltaMag(iC) = Clusters.mMag(iC) - Clusters.aMag(iC);
        Clusters.deltaDays(iC) = days( Catalog.date(idxAftershock) - Catalog.date(idxMainshock) );
%         Clusters.mainVsAftersh_kmDist(iC)   = getDistance( Catalog.lon(idxMainshock), Catalog.lat(idxMainshock), ...
%                                                            Catalog.lon(idxAftershock), Catalog.lat(idxAftershock), 'km' );
        dist_km = getDistance( Catalog.lon(idxMainshock), Catalog.lat(idxMainshock), Catalog.lon(idxAftershock), Catalog.lat(idxAftershock), 'km' );
        Clusters.deltaRupLengths(iC) = dist_km ./ Catalog.rupL_km(min(idxMainshock, idxAftershock));
    
%         if Clusters.deltaRupLengths(iC) > 2.5
%             warning('Test')
%         end
        
    else
        Clusters.triggDepth{iC}     = newCluster.generation;
        Clusters.isLeaf{iC}         = newCluster.isLeaf;
        Clusters.meanLeafDepth(iC)  = mean( newCluster.generation(newCluster.isLeaf) );
        Clusters.maxLeafDepth(iC)   = max( newCluster.generation );
    end
    
    warning on

end