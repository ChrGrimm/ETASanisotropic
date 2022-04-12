function [iC, idxEventsSorted] = merge_clusters( Catalog, Clusters, iC, idxEventsSorted )

    hasSameMainshock = ismember(Clusters.mId, Catalog.id(idxEventsSorted(1)));
        
    if any( hasSameMainshock )
        eventIDs        = unique( [ Catalog.id(idxEventsSorted); ...
                                    cat(1, Clusters.eventIDs{hasSameMainshock}) ] );
        iC              = find(hasSameMainshock);
        idxEvents       = find( ismember(Catalog.id, eventIDs) );
        [~, idxSorted]  = sort(Catalog.mag(idxEvents),'descend');
        idxEventsSorted = idxEvents( idxSorted );
    end
    
end