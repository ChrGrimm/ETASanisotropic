function timeGrid = predefine_interEventTimeGrid( minTimeStep_sec, ...
                                                     maxTimeStep_min, ...
                                                     maxTimeInterval, ...
                                                     c, p)

    % Conversion of time step variables to decimals of a day
    secPerDay           = 24*60*60;
    minPerDay           = 24*60;
    minTimeStep_day     = minTimeStep_sec / secPerDay;
    maxTimeStep_day     = maxTimeStep_min / minPerDay;
    
    % Function handle of Omori law derivative at time t minus at time 0
    funcOmori           = @(t) -p*((t+c).^(-p-1)-c^(-p-1));
    
    % Compue time grid
    minStep_omori       = funcOmori( minTimeStep_day );
    maxStep_omori       = funcOmori( maxTimeInterval );
    omoriGrid           = minStep_omori : minStep_omori : maxStep_omori;
    timeGrid            = nthroot(-omoriGrid/p + c^(-p-1), -p-1) - c;
    timeGrid_diff       = diff(timeGrid);
    if any( timeGrid_diff > maxTimeStep_day )
        idxTooLargeStep = find( timeGrid_diff > maxTimeStep_day, 1, 'first' );
    else
        idxTooLargeStep = length(timeGrid_diff)+1;
    end
    timeGrid            = [ timeGrid(1:idxTooLargeStep), timeGrid(idxTooLargeStep)+maxTimeStep_day : maxTimeStep_day : maxTimeInterval ];
%     timeGrid_sec        = timeGrid * secPerDay;
%     timeGrid_total      = [];
    
end