function [xNewEvents, yNewEvents] = sample_aftershockLocations( BackgrDistr, ...
                                                                TriggerSources, ...
                                                                paramETAS, ...
                                                                nNewEvents, ...
                                                                ModelFuncs, ...
                                                                version )

    %% Initialize output vectors
%     if strcmp(SpatKernel.type, 'none')
    xNewEvents  = zeros(nNewEvents,1);
    yNewEvents  = zeros(nNewEvents,1);
    dist_km     = zeros(nNewEvents,1);
      
    %% Location samples of background events
    if strcmp(version, 'background')
        % Compute cumulative sum of all rates to sample from them
        idxGridP    = randsample( 1:length(BackgrDistr.gridX), nNewEvents, true, BackgrDistr.prob );    
        xNewEvents  = BackgrDistr.gridX(idxGridP) + BackgrDistr.deltaX * (rand(nNewEvents,1)-0.5);
        yNewEvents  = BackgrDistr.gridY(idxGridP) + BackgrDistr.deltaY * (rand(nNewEvents,1)-0.5);
    
    %% Location samples of triggered events
    elseif strcmp(version, 'triggered')
        % Extract data from triggering events
        mag         = TriggerSources.mag;
        xSource     = TriggerSources.x;
        ySource     = TriggerSources.y;
        rupExtent   = TriggerSources.rupExtent;
        strike      = TriggerSources.strike;
        epiPos      = TriggerSources.epiPos;
        typeKernel  = TriggerSources.typeKernel;

        % Sample distance of triggered events to epicenter (isotropic) or rupture segment (anisotropic)
        sigma           = ModelFuncs.fSpatK_width; % @(mag) max( D*exp(gamma*(mag-SimulSettings.M_c)), SpatKernel.minKernelWidth );
        spatRestr       = ModelFuncs.fRestrSpace_degrees(mag,typeKernel);
        q               = paramETAS(8);    
        scalingTrigger  = 1 - ModelFuncs.fSpatK_inn(spatRestr, spatRestr.^2, mag, rupExtent).^(1-q);
        u               = rand(nNewEvents, 1);
        r               = ( sqrt( sigma(mag) .* ( nthroot(1-u.*scalingTrigger, 1-q) - 1) + rupExtent.^2/pi ) - rupExtent/sqrt(pi) ) / sqrt(pi);
        dist_km         = r*111;

        if any(r > spatRestr)
            error('Distance larger than spatial extent sampled!')
        end 
        
        % Sample exact location of isotropically triggered events
        isIso = strcmp( typeKernel, 'iso' ); %mag <  SpatKernel.anisoFromMw;
        if any(isIso)
            [ xNewEvents(isIso), ...
              yNewEvents(isIso) ] = sample_spaceIso( xSource(isIso), ...
                                                     ySource(isIso), ...
                                                     r(isIso) );
        end
        
        % Sample exact location of anisotropically triggered events
        isAniso = strcmp( typeKernel, 'aniso' ); %mag >= SpatKernel.anisoFromMw;
        if any(isAniso)
            [ xNewEvents(isAniso), ...
              yNewEvents(isAniso) ] = sample_spaceAniso( xSource(isAniso), ...
                                                         ySource(isAniso), ...
                                                         rupExtent(isAniso), ...
                                                         strike(isAniso), ...
                                                         epiPos(isAniso), ...
                                                         r(isAniso) );
        end
        
    end
    
end
