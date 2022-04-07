function [ TargetSettings, ...
            ModelFuncs ] = preprocess_targetSettings( TargetSettings, ...
                                                      ModelFuncs )

    %% Target polygon
    % Load and store lon-lat polygons
    load(TargetSettings.pathPolygon, 'polygonTarget', 'polygonComple')
    TargetSettings.polygonTarget = check_polygonInput( polygonTarget ); 
    TargetSettings.polygonComple = check_polygonInput( polygonComple );
    
    % Determine conversion formulas for lon-lat-coordinates to x-y space
    ModelFuncs = compute_conversionsLonLat2XY( ModelFuncs, polygonTarget );
    
    % Convert target polygon
    TargetSettings.polygonXY(:,1) = ModelFuncs.fConvert_lon2x( TargetSettings.polygonTarget(:,1) );
    TargetSettings.polygonXY(:,2) = ModelFuncs.fConvert_lat2y( TargetSettings.polygonTarget(:,2) );
%     TargetSettings.polygonCenter  = [lonCenter, latCenter];
    
    % Compute area inside of polygon
    TargetSettings.polyArea      = polyarea( TargetSettings.polygonXY(:,1), TargetSettings.polygonXY(:,2) );

    %% ETAS estimation
    if isfield(TargetSettings, 'tStart')
        % Compile time window (dates)
        TargetSettings.tWindow = [ TargetSettings.tStart, TargetSettings.tEnd ];
        % Convert time window to days format
        TargetSettings.tWindow_days  = days( TargetSettings.tWindow - TargetSettings.tWindow(1,1) ); %[ 0, days(diff(TargetSettings.timeWindow)) ];
        
        if any( TargetSettings.tWindow(:,2) < TargetSettings.tWindow(:,1) )
            error('Start and end date of target period are in wrong chronological order')
        end
        % History start must be before start of target time window
        if TargetSettings.tWindow(1,1) < TargetSettings.tIni 
            error('Initialization date is after target window start date.')
        end   
        
    %% ETAS simulation
    elseif isfield(TargetSettings, 'nYearsPerPeriod')
        TargetSettings.tWindow_days = [0, 365*TargetSettings.nYearsPerPeriod];
    else
        error('Missing inputs')
    end
    
end
