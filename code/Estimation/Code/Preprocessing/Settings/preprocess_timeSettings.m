function [TimeSettings, ModelFuncs] = preprocess_timeSettings( TimeSettings, ...
                                                                TargetSettings, ...
                                                                ModelFuncs )

    %% Modify time settings
    % Cap temporal restriction at length of time window
    TimeSettings.restr_days = min( TimeSettings.restr_days, TargetSettings.tWindow_days(end,2) );
    
    %% Temporal restriction
    ModelFuncs.fRestrTime_days = @(mag) TimeSettings.restr_days * ones(size(mag));
    
end