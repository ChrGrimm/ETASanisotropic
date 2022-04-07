function [ completeTimeWindows, ...
           isInCompleteTime ] = apply_shortTermIncompl2timeWindow( Catalog, ...
                                                                    TargetSettings, ...
                                                                    ModelSettings )

%     isInCompleteTime = false(size(Catalog.t)); 
    Catalog.t       = days( Catalog.date-TargetSettings.tWindow(1,1) );
    tIncomplPeriods = zeros(length(Catalog.mag),2);
    nRows           = 1;
    Mc              = TargetSettings.Mc;
    % isActive = false(length(Catalog.mag),1);
    
    tWindow_days    = TargetSettings.tWindow_days;
    idxStaiEvents   = find( Catalog.mag >= ModelSettings.staiFromMw & ...
                            Catalog.t >= tWindow_days(1,1) & Catalog.t <= tWindow_days(end,2) & ...
                            inpolygon(Catalog.lon, Catalog.lat, TargetSettings.polygonTarget(:,1), TargetSettings.polygonTarget(:,2)) )';

    for iEv=idxStaiEvents

        %% Compute start and end time of incompleteness period
        % Start time is the event occurrence time
        tStartIncompl  = Catalog.t(iEv);
        % End time is dependent on chosen incompleteness function
        switch ModelSettings.typeSTAI
            case 'Helmstetter2006'
                tEndIncompl    = Catalog.t(iEv) + 10.^((Catalog.mag(iEv)-4.5-Mc)/0.75);
            case 'RidgecrestMw64'
                tEndIncompl    = Catalog.t(iEv) + 10.^(Catalog.mag(iEv)-5-Mc);
            otherwise
                error('Unknown type of short-term aftershock incompleteness!')
        end

        %% Update matrix storing incomplete periods
        if tEndIncompl <= tWindow_days(1)
            % Current incompleteness period ends before the start of the target time window
            continue
        elseif any( tEndIncompl <= tIncomplPeriods(1:nRows-1,2) ) % tStartIncompl < tIncomplPeriods(1:nRows-1,2) &
            % Current incompleteness period is embedded in another (that starts earlier and ends later)
            continue
        elseif any( tStartIncompl <= tIncomplPeriods(1:nRows-1,2) & tEndIncompl > tIncomplPeriods(1:nRows-1,2) )
            % Current incompleteness period starts before, but ends after the end of another
            % The two incompleteness periods are merged (by updating the end of the period)
            thisIncomplPeriod                       = tStartIncompl < tIncomplPeriods(1:nRows-1,2) & ...
                                                      tEndIncompl   > tIncomplPeriods(1:nRows-1,2);
            tIncomplPeriods( thisIncomplPeriod, 2 ) = tEndIncompl;
        elseif all( tStartIncompl > tIncomplPeriods(1:nRows-1,2) )
            % Current incompleteness period starts after the end of any others
            % New incompleteness period is added to dataset
            tIncomplPeriods(nRows,:)    = [tStartIncompl, tEndIncompl];
            nRows                       = nRows+1;
        end

    end

    % Delete rows of zeros
    tIncomplPeriods(nRows:end,:) = [];
    if isempty(tIncomplPeriods)
        completeTimeWindows = tWindow_days;
    else
        completeTimeWindows = [ tWindow_days(1,1), tIncomplPeriods(1,1); ...
                                tIncomplPeriods(1:end-1,2), tIncomplPeriods(2:end,1); ...
                                tIncomplPeriods(end,2), tWindow_days(end,2) ];
    end
    isInCompleteTime = any( Catalog.t >= completeTimeWindows(:,1)' & Catalog.t <= completeTimeWindows(:,2)', 2 );

end