function isActiveFault = evaluate_activeFaultsPerPeriod( MET, Faults, EventListWPET )

    %% Preallocate output boolean matrix
    nPeriods        = max(EventListWPET.periodID);
    nFaultGridP     = size(Faults, 1);
    isActiveFault   = true(nFaultGridP, nPeriods);
    
    %% Set to false, if characteristic fault already ruptured in period
    for iPeriod = 1:nPeriods
        % Determine already ruptured characteristic faults
        eventIDs    = EventListWPET.id(EventListWPET.periodID==iPeriod);
        isCharEvent = ismember(MET.eventID, eventIDs) & MET.isCharEvent;
        charFaultID = MET.faultID(isCharEvent);
        
        % Set already ruptured characteristic faults to "false" (cannot rupture again)
        isActiveFault(ismember(Faults.faultID, charFaultID), iPeriod) = false;
    end
    
end