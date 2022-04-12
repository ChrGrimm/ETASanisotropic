function TriggerRelations = initialize_tTriggerRelations( Catalog, strModel )

    
    %% Extract target events (for ETAS version)
    if strcmp( strModel, 'ETAS' )
        isTarget    = Catalog.flag > 0 & ~Catalog.isDupl;
        Catalog     = Catalog( isTarget, : );
    end
    
    %% Initialize table
    nEvents                     = size(Catalog,1);
    TriggerRelations            = table;  
    TriggerRelations.id         = Catalog.id;
    TriggerRelations.parentID   = -1*ones(nEvents,1);
    
    if strcmp( strModel, 'ETAS' )
        TriggerRelations.backgrProb         = Catalog.backgrProb;
        TriggerRelations.triggProb          = -1*ones(nEvents,1);
        
    elseif strcmp( strModel, 'ZBZ' )
        TriggerRelations.distance_ZBZ   = zeros(nEvents,1);
        TriggerRelations.tij_rescaled   = zeros(nEvents,1);
        TriggerRelations.rij_rescaled   = zeros(nEvents,1);   
        
    end
    
%     TriggerRelations.offsprings     = cell(nEvents,1);
    TriggerRelations.nOffspr_count  = zeros(nEvents,1);

end