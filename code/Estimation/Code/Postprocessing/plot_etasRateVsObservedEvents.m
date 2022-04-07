function plot_etasRateVsObservedEvents( Inputs, Catalog, paramETAS, SpatData, TempData, ModelFuncs, nBackgrPerDay, pathResultsFolder )
% This function plots the integrated ETAS event rate against the observed number of events in 7-days or 
% hourly time intervals (depending on catalog length)
% Note: This plotting function can only be called from inside of the ETAS estimation, as it requires 
% large size data inputs that are not stored in the results

    %% Preprocessing
    % Load estimation results
    [mu,~,~,~,p,D,gamma,q,Tb]   = extract_paramETASvalues( paramETAS, 'default' );
    fTempK_inn                  = ModelFuncs.fTempK_inn;
    
    % Define time grid
    if diff(Inputs.TargetSettings.tWindow_days(1,:)) > 100
        stepSize_days = 7;
    else
        stepSize_days = 1/24;
    end
    timeGrid        = Inputs.TargetSettings.tWindow_days(1,1):stepSize_days:Inputs.TargetSettings.tWindow_days(end,2);
    nTimeGridWindow = length(timeGrid)-1;
    
    %% Evaluate integrated event rates between time grid points
    % Preallocate rate vectors
    triggerRate = zeros( nTimeGridWindow, 1 );
    backgrRate  = zeros( nTimeGridWindow, 1 );
    
    %% Precomputations
    % Aftershock Productivity
    Productivity = ModelFuncs.fProductivity(Catalog.mag, Inputs.TargetSettings.Mc);
    % Spatial integral
    SpatIntegr = estimate_spatialIntegral( 0, Catalog, Inputs, [D, gamma, q], ModelFuncs, SpatData );
    % for ETAS-Incomplete only
    if strcmp(Inputs.ModelSettings.modelType, 'ETAS-Incomplete')
        N0              = estimate_N0( Catalog, Inputs, TempData, nBackgrPerDay/mu, sqrt(paramETAS), ModelFuncs, SpatIntegr, -1, -1, -1, 'no gradient' );
        expMinusN0      = exp(-N0)';
        timeGrid_etasi  = TempData.timeGrid;
    end
    
    % Loop over each time interval
    for iWindow = 1:nTimeGridWindow
        % Evaluate start and end point of current time grid window
        tStart = timeGrid(iWindow);
        tEnd   = timeGrid(iWindow+1);
        
        % Estimate background rate
        backgrRate(iWindow)  = sum(tEnd-tStart) * nBackgrPerDay; 
        
        % Estimate trigger rate (depending on model type
        if strcmp(Inputs.ModelSettings.modelType, 'ETAS')
            TempIntegr          = estimate_temporalIntegral( 0, Catalog, p, fTempK_inn, 1:N, [tStart,tEnd] );
            triggerRate(iWindow)= sum(Productivity .* SpatIntegr .* TempIntegr);
            
        elseif strcmp(Inputs.ModelSettings.modelType, 'ETAS-Incomplete')
            idxThisWindow       = timeGrid_etasi>tStart & timeGrid_etasi<tEnd;
            thisGrid            = [tStart; timeGrid_etasi(idxThisWindow); tEnd];
            thisExpMinusN0      = [ interp1(timeGrid_etasi, expMinusN0, tStart, 'linear', 'extrap'); ...
                                    expMinusN0(idxThisWindow); ...
                                    interp1(timeGrid_etasi, expMinusN0, tEnd, 'linear', 'extrap') ];
            triggerRate(iWindow)= (tEnd-tStart)/Tb - trapz( thisGrid, thisExpMinusN0 ) / Tb - backgrRate(iWindow);
        end
        
    end
    
    %% Sum up triggered and background rates
    rate_integrated = triggerRate + backgrRate;
    coeffOfVar_etas = std(rate_integrated)/mean(rate_integrated);

    %% Count observed event rates
    rate_observed   = histcounts( Catalog.t(~Catalog.isDupl), timeGrid );
    coeffOfVar_obs  = std(rate_observed)/mean(rate_observed);
        
    %% Evaluate number of observed events between time grid points
    timeGridPoints  = (timeGrid(2:end)+timeGrid(1:end-1))/2;
    plot(timeGridPoints, rate_observed, 'k-', 'LineWidth', 1.5)
    hold on
    plot(timeGridPoints, rate_integrated, 'r-', 'LineWidth', 1.5)
    h_leg = legend('observed event rate', 'integrated ETAS rate');
    set(h_leg,'Location','NorthEast','Interpreter','none')
    
    if stepSize_days == 7
        strTitle_line1 = 'ETAS vs. Observed Event Rates (per 7 days)';
        
    elseif stepSize_days == 1/24
        strTitle_line1 = 'ETAS vs. Observed Event Rates (per hour)';
        
    else
        warning('Unknown step size')
    end
    
    strTitle_line2 = ['Coefficient of Variation: ', num2str(coeffOfVar_etas), ' (estimated) vs. ', num2str(coeffOfVar_obs), ' (observed)'];
    title( {['\fontsize{12}', strTitle_line1], ['\fontsize{9}', strTitle_line2]} )
    xlabel('Time t (in days since model start date)')
    ylabel('Number of events per time interval')
    
end    