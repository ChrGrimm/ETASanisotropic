function check_clusters( Clusters, TriggerRel, Catalog )

    allEventIDs = [];
    for i=1:size(Clusters,1)
        allEventIDs = [allEventIDs; Clusters.eventIDs{i}];
    end
    
    if ~all( sort(allEventIDs)==Catalog.id )
        error('Every event ID must appear exactly ones in Clusters!')
    end
    
    if ~all( TriggerRel.id == Catalog.id( Catalog.flag>0 ) )
        error('TriggerRelationsETAS should contain all target events!')
    end

end