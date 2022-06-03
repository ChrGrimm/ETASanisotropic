function isTriggered = traceBack_branchingTree( triggerID, id, iEv, refEv )

    %% Initializations
    % Output boolean vector
    isTriggered = false(size(refEv));
    % ID of event under observation
    iID         = id(iEv);
    
    %% Trace back branching tree for each reference event
    % Loop over reference events
    for i=1:length(refEv)
        % current reference event
        iRefEv          = refEv(i);
        % ID of parent event and IDs of already investigated reference events
        nextTriggerID   = triggerID(iRefEv);
        oldRefIDs       = id(refEv(1:(i-1)));
        
        %% Check if events have triggering relation
        if ismember(nextTriggerID, oldRefIDs)
            % If current reference event triggered by an already investigated reference event...
            isTriggered(i) = isTriggered(nextTriggerID==oldRefIDs);
%         elseif nextTriggerID == 0 || nextTriggerID == -1
%             % If reference event is background or complementary
%             isTriggered(i) = false;
        else
            %% Trace back branching tree
            while ~isempty(nextTriggerID)
                if nextTriggerID == iID
                    % If trigger event equals event under consideration
                    isTriggered(i) = true;
                    break;            
                else 
                    % Update trigger ID
                    nextTriggerID = triggerID( id==nextTriggerID );
                end
            end
        end
    end
    
end