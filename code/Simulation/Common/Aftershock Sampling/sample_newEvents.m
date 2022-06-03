function [ SynthEventSet, ...
            idxTriggerEvents ] = sample_newEvents( SynthEventSet, ...
                                                   BackgrDistr, ...
                                                   ExperimentSettings, ...
                                                   TargetSettings, ...
                                                   SpaceSettings, ...
                                                   TriggerSources, ...
                                                   paramETAS, ...
                                                   ModelFuncs, ...
                                                   InputsWPET, ...
                                                   version )
 
    %% Compute number of new events
    if strcmp(version, 'background')
        % Sample number of background events
        nNewEvents = poissrnd( BackgrDistr.nBackgrMean );
    elseif strcmp(version, 'triggered')
        % Number of aftershocks = Number of rows of TriggerSources
        nNewEvents = size(TriggerSources, 1);
    end
    
    %% Sample occurrence times of new events
    [t, TriggerSources] = sample_aftershockTimes( TriggerSources, ...
                                                    TargetSettings.tWindow_days, ...
                                                    paramETAS, ...
                                                    ModelFuncs.fTempK_inn, ...
                                                    ModelFuncs.fRestrTime_days, ...
                                                    nNewEvents, ...
                                                    version);

    %% Stop function if no new events are remaining
    nNewEvents = length(t);
    if nNewEvents == 0 
        idxTriggerEvents = [];
        return;
    end
    
    %% Sample new event locations
    [x, y] = sample_aftershockLocations( BackgrDistr, ...
                                         TriggerSources, ...
                                         paramETAS, ...
                                         nNewEvents, ...
                                         ModelFuncs, ...
                                         version );
                                             
    isInPoly = inpolygon(x, y, TargetSettings.polygonXY(:,1), TargetSettings.polygonXY(:,2));
    x(~isInPoly) = [];
    y(~isInPoly) = [];
    t(~isInPoly) = [];
    TriggerSources(~isInPoly,:) = [];
    
    nNewEvents = length(x);
    if nNewEvents == 0 
        idxTriggerEvents = [];
        return;
    end    
    
    %% Sample aftershock magnitudes
    if strcmp(ExperimentSettings.TypeSimulation, 'Catalog')
        mag = sample_magnitudesByGR( nNewEvents, ...
                                     ExperimentSettings, ...
                                     TargetSettings.Mc, ...
                                     TargetSettings.Mmax, ...
                                     paramETAS(10) );   
                                 
    elseif strcmp(ExperimentSettings.TypeSimulation, 'WPET')  
        
        if isempty(ExperimentSettings.sampleMagnitudes) || isnumeric(ExperimentSettings.sampleMagnitudes)
            mag = sample_magnitudesByGR( nNewEvents, ...
                                         ExperimentSettings, ...
                                         TargetSettings.Mc, ...
                                         TriggerSources.mag, ...
                                         paramETAS(10) );   

        elseif strcmp(ExperimentSettings.sampleMagnitudes, 'MET')
            mag = sample_magnitudesByMET( x, y, InputsWPET, TriggerSources.periodID );

        end
        
    end
    
    %% Determine and sample spatial event characteristics
    % Determine type of kernel
    isAniso             = mag >= SpaceSettings.anisoFromMw;
    typeKernel          = cell(nNewEvents,1);
    typeKernel(:)       = {'iso'};
    typeKernel(isAniso) = {'aniso'};
    nAniso              = sum(isAniso);
    
    % Determine rupture extent and spatial restriction
    rupExtent           = zeros(nNewEvents,1);
    rupExtent(isAniso)  = ModelFuncs.fRupLength( mag(isAniso) );
    spatRestr           = ModelFuncs.fRestrSpace_degrees(mag,typeKernel);

    %% Match event ID from MET
    if strcmp(ExperimentSettings.TypeSimulation, 'WPET')    
        [id, depth, strike, epiPos] = match_aftershocks2met( x, y, mag, InputsWPET.MET );
    else
        id                  = -1 * ones(nNewEvents,1);
        depth               = -1 * ones(nNewEvents,1);     
        % Sample strikes and epicenter location
        strike              = nan(nNewEvents,1);
        strike(isAniso)     = 360*rand(nAniso,1);
        epiPos              = nan(nNewEvents,1);
        epiPos(isAniso)     = rand(nAniso,1);
    end
            
    %% Store trigger and cluster information
    if strcmp(version, 'background')
        % New events are not triggered and start their own, new cluster
        periodID    = -1*ones(nNewEvents,1);
        triggerID   = zeros(nNewEvents,1);
        clusterID   = size(SynthEventSet,1) + (1:nNewEvents)';
    elseif strcmp(version, 'triggered')
        % New events are triggered and belong to the cluster of the triggering event
        periodID    = TriggerSources.periodID;
        triggerID   = TriggerSources.id;
        clusterID   = TriggerSources.clusterID;
    end
    
    %% Fill other columns of new events
    % Event ID
%     id              = SynthEventSet.id(end) + (1:nNewEvents)';
    % Event flag
    if strcmp(version, 'background')
        flag            = ones(nNewEvents,1);
    elseif strcmp(version, 'triggered')
        [inPol, onPol]  = inpolygon(x, y, TargetSettings.polygonXY(:,1), TargetSettings.polygonXY(:,2));
        flag            = 2*inPol; % - 0.5*onPol;
    end
    % Longitude, Latitude, Depth, Event Weight
    evWeight        = ones(nNewEvents,1);
    idxMainsh       = TriggerSources.idxMainsh;

    %% Add new events to synthetic event set table
    newData             = table(periodID, id, t, x, y, depth, mag, flag, evWeight, typeKernel, rupExtent, spatRestr, ...
                                strike, epiPos, triggerID, clusterID, idxMainsh);
    idxTriggerEvents    = size(SynthEventSet,1) + (1:nNewEvents)'; 
    SynthEventSet       = [SynthEventSet; newData];
    
end