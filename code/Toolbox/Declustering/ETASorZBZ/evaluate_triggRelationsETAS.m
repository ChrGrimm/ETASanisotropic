function TriggerRelations = evaluate_triggRelationsETAS( Catalog, triggProbMatrix )

    %% Compile table "TriggerRelations"
    TriggerRelations = initialize_tTriggerRelations( Catalog, 'ETAS' );
    
    %% Evaluate parent event, if triggered
    nEvents                     = size(TriggerRelations,1);
    
    for iEv = 1:nEvents
        if TriggerRelations.backgrProb(iEv) >= 0.5
            continue;
        else
            [ ~, idxParent ]                = max(triggProbMatrix(:,iEv));
            TriggerRelations.parentID(iEv)  = Catalog.id(idxParent');
            TriggerRelations.triggProb(iEv) = triggProbMatrix(idxParent,iEv);
        end
    end
    
    % Evaluate offsprings
    for iEv = 1:nEvents
        TriggerRelations.nOffspr_count(iEv) = sum( ismember(TriggerRelations.parentID, TriggerRelations.id(iEv)) );
    end
    isTarget                            = Catalog.flag > 0;
    TriggerRelations.nOffspr_sumProb    = sum(triggProbMatrix(isTarget,:), 2);
      
end