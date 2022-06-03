function [ doubletDegree, ...
           doubletType, ...
           extraRowsDoubletType ] = evaluate_doubletDegreeAndType( idxTimeSpace, ...
                                                                   iEvent, ...
                                                                   mag, ...
                                                                   magnWindow, ...
                                                                   flag, ...
                                                                   triggerID, ...
                                                                   eventID, ...
                                                                   clusterID )
         
    %% Evaluate doublet degree
    % Apply magnitude window for doublets
    isInMagn        = round( abs(mag(idxTimeSpace)-mag(iEvent)), 1 ) <= magnWindow;

    % Evaluate number of doublet partners for current event
    doubletDegree   = sum(isInMagn);

    % If doublet given, evaluate relation of doublet partners
    if doubletDegree && flag(iEvent)>-1                                                         

        refEv           = idxTimeSpace( isInMagn );
        isTriggered     = traceBack_branchingTree( triggerID, eventID, iEvent, refEv );
        relationScore   = sum( [triggerID(refEv)==eventID(iEvent), ...
                                isTriggered, ...
                                clusterID(refEv)==clusterID(iEvent)], 2 );
    %     relationScore   = sum([triggerID(refEv)==eventID(iEvent), clusterID(refEv)==clusterID(iEvent)],2);
    
        if doubletDegree == 1
            
            doubletType          = [mag(iEvent), relationScore];
            extraRowsDoubletType = NaN(0,2);
            
        else
            
            doubletType          = [ mag(iEvent), relationScore(1) ];
            nAdd                 = length(relationScore)-1;
            extraRowsDoubletType = [ mag(iEvent)*ones(nAdd,1), relationScore(2:end) ];
            
        end
        
    else
        
        doubletType             = NaN(1,2);
        extraRowsDoubletType    = NaN(0,2);
        
    end

end