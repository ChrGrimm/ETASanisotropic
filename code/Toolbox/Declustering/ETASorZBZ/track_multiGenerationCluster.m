function Cluster = track_aftershockCascade( TriggerRelations, Catalog, iEvent )

    warning off
    
    %% Preallocate new cluster
    Cluster             = table;
    Cluster.evIDs       = Catalog.id(iEvent);
    Cluster.generation  = 0;
    Cluster.isLeaf      = false;

    %% Preallocate auxiliary variables
    % Trigger generation counter
    iGeneration         = 0;
    % Cluster size
    clusterSize         = 1;
    % Cluster primarily has one source
    triggeringIDs       = Catalog.id(iEvent);

    %% Loop over trigger generations
    while true
        % Update trigger generation counter
        iGeneration         = iGeneration + 1;

        % Find new generation of triggered events
        offspringIDs    = TriggerRelations.id( ismember(TriggerRelations.parentID, triggeringIDs) );

        % Leave loop of no triggered events
        if isempty( offspringIDs )
            break;
        end
        
        % Row indices to fill with new trigger generation
        idxLow              = clusterSize + 1;
        idxUpp              = clusterSize + length(offspringIDs);
        clusterSize         = idxUpp; 

        % Update vector of cluster events
        Cluster.evIDs(idxLow:idxUpp)        = offspringIDs;
        % Update vector of trigger generation markers
        Cluster.generation(idxLow:idxUpp)   = iGeneration;
        
        % Update vector of leaf markers
        Cluster.isLeaf(idxLow:idxUpp)       = ~ismember(offspringIDs, TriggerRelations.parentID);  

        % Update vector of last generation event IDs
        triggeringIDs   = offspringIDs;

    end
    
    warning on

end