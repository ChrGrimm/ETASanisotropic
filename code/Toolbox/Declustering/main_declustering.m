function [ Clusters, ...
           Catalog ] = main_declustering( Catalog, ...
                                            TriggerRelations, ...
                                            Settings, ...
                                            declusterMethod )

    if ~ismember(declusterMethod, {'window method', 'ETAS', 'ZBZ'})
        error('Unknown method: Choose from ''window method'', ''ETAS'', ''ZBZ'':')
    end
    
    %% For windowing method: Sort catalog by descending magnitudes
    if strcmp(declusterMethod, 'window method')
        Catalog             = filter_origCatalog( Catalog, Settings.CatalogSettings, 'declustering' );
        Catalog             = sortrows(Catalog, 'mag', 'descend');
        Catalog.rupL_km     = estimate_rupSize( Catalog.mag, ...
                                                Settings.WindowSettings.tecType, ....
                                                Settings.WindowSettings.faultStyle, ...
                                                'km' );
        Catalog.spatRestr   = Settings.WindowSettings.xTimesRupL * Catalog.rupL_km;
    end
    
    %% Identify clusters
    % Initialize variables
    Clusters        = table;
    nEvents         = size(Catalog, 1);
    isPartOfCluster = false(nEvents, 1);
    
    % Loop over all events
    for iEvent = 1:nEvents

        if isPartOfCluster(iEvent)
            % Skip events that are already assigned to another cluster
            continue;
            
        else
            
            %% Track cluster initiated by current event
            if strcmp(declusterMethod, 'window method')
                newCluster = identify_eventsInSpaceTimeWindow( Settings.WindowSettings, Catalog, iEvent );
            else
                newCluster = track_multiGenerationCluster( TriggerRelations, Catalog, iEvent ); 
            end

            if size(newCluster,1) > 1
                % Mark events as part of cluster
                isPartOfCluster( ismember(Catalog.id, newCluster.evIDs) ) = true;
                
                % Update cluster table
                Clusters = update_clusterStatistics( Clusters, newCluster, Catalog, declusterMethod );   
            end

        end
    end
    
    %% Add Singles
    nSingles                        = sum(~isPartOfCluster);
    idxSingles                      = size(Clusters,1) + (1:nSingles)';
    warning off
    Clusters.clusterID(idxSingles)  = idxSingles;
    Clusters.size(idxSingles)       = 1;
    Clusters.mId(idxSingles)        = Catalog.id(~isPartOfCluster);
    Clusters.mMag(idxSingles)       = Catalog.mag(~isPartOfCluster);
    Clusters.mDate(idxSingles)      = Catalog.date(~isPartOfCluster);
    warning on
    
    %% Order clusters by magnitude
    Clusters                        = sortrows(Clusters, {'size', 'mMag'}, 'descend');
    Clusters.clusterID              = (1:size(Clusters,1))';
    
    %% Plausibility checks on results
%     check_clusters( Clusters, TriggerRelations, Catalog )
    
end