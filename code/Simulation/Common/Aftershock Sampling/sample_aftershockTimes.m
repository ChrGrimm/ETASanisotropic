function [t, TriggerSources] = sample_aftershockTimes( TriggerSources, ...
                                                         tWindow_days, ...
                                                         paramETAS, ...
                                                         fTempKernel_inn, ...
                                                         fRestrTime_days, ...
                                                         nNewEvents, ...
                                                         version)
                                 
    if strcmp(version, 'background')
        % Uniformly distributed over simulated time window
        t       = tWindow_days(2) * rand(nNewEvents, 1);
        
    elseif strcmp(version, 'triggered')
        % Distributed according to normalized Omori-Utsu Law, solved for t
        c           = paramETAS(4);
        p           = paramETAS(5);
        T           = fRestrTime_days( TriggerSources.mag );
        tempScaling = 1/(1-p) * ( fTempKernel_inn(T).^(1-p) - fTempKernel_inn(0).^(1-p) );
        diffDays    = nthroot( (1-p)*rand(nNewEvents,1).*tempScaling + c^(1-p), 1-p ) - c;
        t           = TriggerSources.t + diffDays;
        
        % Remove events that occurred outside of time window
        isInTimeWindow                      = t >= 0 & t <= tWindow_days(2);
        t(~isInTimeWindow)                  = [];
        diffDays(~isInTimeWindow)           = [];
        TriggerSources(~isInTimeWindow, :)  = [];
        
    end
    
end