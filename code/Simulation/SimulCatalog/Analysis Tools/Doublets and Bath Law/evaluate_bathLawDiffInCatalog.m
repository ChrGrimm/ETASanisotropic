function [ bathLawResults, ...
            eventRelation ] = evaluate_bathLawDiffInCatalog( Catalog, ...
                                                             timeWindow, ...
                                                             spaceUnit, ...
                                                             isSynthCat )
    %% FUNCTION TO DELETE
    % Store often used columns
    x           = Catalog.x;
    y           = Catalog.y;
    t           = Catalog.t;
    mag         = Catalog.mag;
    rupL        = Catalog.rupL;
    id          = Catalog.id;
    if isSynthCat
        triggerID   = Catalog.triggerID;
        clusterID   = Catalog.clusterID;
    end
    % Extract event indices for which bath law should be evaluated
    idxEvents   = find(Catalog.mag>=5.5 & abs(Catalog.flag)>0.1 & Catalog.t <= timeWindow(end,2)-365);
    % Initialize output and auxiliary variable
    pairs       = zeros(length(idxEvents),2);
    eventRelation = zeros(length(idxEvents),2);
    toDelete    = [];
    
    isInTime    = t-t(idxEvents)'>0 & t-t(idxEvents)'<=365;

    %% Find smallest magnitude difference
    for k = 1:length(idxEvents)
        % Consider event iEv
        iEv         = idxEvents(k);
        % Find all events within the time window
        idxInTime   = find(isInTime(:,k));
        % Find all events within the space window
        isInSpace   = getDistance( x(idxInTime), y(idxInTime), x(iEv), y(iEv), spaceUnit) <= 2.5*rupL(iEv);
        if any(isInSpace)
            % Find maximum aftershock magnitude and associated event ID
            idxInTimeAndSpace       = idxInTime(isInSpace);
            [maxAftershock, idxMax] = max( mag(idxInTimeAndSpace) );
            % Store magnitude difference
            pairs(k,:)              = [mag(iEv), mag(iEv)-maxAftershock];            
            if isSynthCat
                % Evaluate relation of the two events
                refEv                   = idxInTimeAndSpace(idxMax);
                relationScore           = sum([triggerID(refEv)==id(iEv), clusterID(refEv)==clusterID(iEv)]);
                eventRelation(k,:)      = [mag(iEv), relationScore];
            end
        else
            toDelete = [toDelete; k];
%             pairs(k,:)  = [mag(iEv), mag(iEv)-min(mag)+0.1];
        end

    end

    pairs(toDelete,:)           = [];
    eventRelation(toDelete,:)   = [];

    [ bathLawResults, ~, idx ]  = unique(pairs, 'rows');
    bathLawResults(:,3)         = sum(unique(idx)==idx', 2);

end
        