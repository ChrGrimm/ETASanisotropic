function mag = sample_magnitudesByMET( xSampled, ySampled, MET, Zones ) %, Faults, periodID )

    %% Preallocation
    mag = zeros(length(xSampled), 1);
    
    %% Compute distances to closest 
    zoneLocations   = table2array(Zones(:,{'x','y'}));
%     faultLocations  = table2array(Faults(:,{'x','y'}));

    [idxZoneGridP, dist2zone]   = knnsearch(zoneLocations, [xSampled,ySampled]);
%     [~, dist2fault]             = knnsearch(faultLocations, [xSampled,ySampled]);
    
    %% Sample zone event magnitudes
    probZoneGridP   = zeros(max(idxZoneGridP),1);
    for iZoneGridP = unique(idxZoneGridP)'
        idxMET              = find(MET.idZoneGridP == iZoneGridP);
        idxAftershocks      = find( idxZoneGridP==iZoneGridP );
        mag(idxAftershocks) = randsample( MET.mag(idxMET), length(idxAftershocks), true, MET.eventProb(idxMET) );
        
        % Compute annual probability of events at this zone grid point
%         probZoneGridP(iZoneGridP) = sum(MET.eventProb(idxMET));
    end
    
%     %% Sample zone or fault event magnitudes
%     isClosest2fault = dist2fault < dist2zone;
% %     idxZoneEvent    = find(isZoneEvent);
%     idxClosest2fault = find(~isClosest2fault);
%     
%     for iEvent = idxClosest2fault'
%         % Recompute distance to closest active fault (given the period)
%         iPeriod                     = periodID(iEvent);
%         idxActiveFaults             = find( InputsWPET.isActiveFault(:,iPeriod) );
%         [idxClosestFault, dist2fault] = knnsearch(faultLocations(idxActiveFaults,:), [xSampled(iEvent),ySampled(iEvent)]);
%         % If distance to closest, active fault is smaller than distance to nearest one point
%         if dist2fault < dist2zone(iEvent)
%             % Draw MET event on nearest fault, depending on sum of fault event probabilities
%             idxThisFault = idxActiveFaults(idxClosestFault);
%             idxMET          = find(MET.idFaultGridP == idxThisFault);
%             probFaultGridP  = sum(MET.eventProb(idxMET));
%             if rand <= probFaultGridP/(probZoneGridP(iEvent)+probFaultGridP)
%                 mag(iEvent) = randsample( MET.mag(idxMET), 1, true, MET.eventProb(idxMET) );
%                 InputsWPET.isActiveFault(idxThisFault,iPeriod) = false;
%             end
%         end
%         
%     end

end