function eventLoss = sample_lossesFromELT( WPET_table, ELT )

    %% Extract columns from tables
    sampleEventIDs  = WPET_table.id;
    eltEventIDs     = ELT.eventID_agg;
    eltEventProbs   = ELT.eventProb;
    eltLoss         = ELT.loss;
    
    %% Determine unique sampled event IDs
    uniqueIDs       = unique(sampleEventIDs);
    
    %% Match losses from ELT to sampled events
    % Initialize output vector with zero losses
    eventLoss = zeros(size(WPET_table, 1), 1);
    
    for iUniID = 1:length(uniqueIDs)
        % Find current event ID in sampled event set and ELT
        thisID      = uniqueIDs(iUniID);
        isSampled   = sampleEventIDs==thisID;
        isELT       = eltEventIDs==thisID;
        
        % Match loss (if any match with ELT - otherwise keep loss=0)
        nMatches = sum(isELT);
        if nMatches==1
            % If exactly one match, take that loss
            eventLoss(isSampled) = eltLoss(isELT);
        elseif nMatches>1
            % If more than one match, sample loss according to event
            % probabilities
            eventLoss(isSampled) = randsample(eltLoss(isELT), sum(isSampled), true, eltEventProbs(isELT));
        end
    end

end