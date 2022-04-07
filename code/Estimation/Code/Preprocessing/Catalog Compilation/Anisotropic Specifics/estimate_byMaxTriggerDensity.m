function [ strikes, ...
           epiPos ] = estimate_byMaxTriggerDensity( sampleStrike, ...
                                                     sampleEpiPos, ...   
                                                     Catalog, ...
                                                     Inputs, ...
                                                     ModelFuncs, ...
                                                     idxAniso )
%
% INPUTS:
% 
% sample:   numerical array     (NxS)   sample of nodal plane solutions or epicenter positions;
%                                       1st row: strike or epicenter position sample
%                                       2nd row: dip sample
%                                       One column pair of strike and dip represents a unique sample                                       
   
    %% Define settings
    % Time window in which events are counted to maximize trigger rate  
    
    %% Preprocess input data
    % Extract parameters from vector paramETAS
    [~,~,~,~,p,~,~,q] = extract_paramETASvalues( Inputs.ModelSettings.iniParamETAS, 'default' );
    
    % Extract space unit
    spaceUnit   = 'degree'; %Inputs.SpaceSettings.spaceUnit;
    minDist     = Inputs.SpaceSettings.minXeventDist;
    
    % Extract event information from columns of table Catalog
    x           = Catalog.x;
    y           = Catalog.y;                                          
    t           = Catalog.t;
    mag         = Catalog.mag;
    spatRestr   = Catalog.spatRestr;  
    rupExtent   = Catalog.rupExtent;
    
    %% Initializations of output % auxiliary variables
    % Matrix of distributed trigger rates per sample item:
    % 1st dimension: number of anisotropic events
    % 2nd dimension: size of strike/epicenter location sample
    nEvents         = size( sampleStrike, 1 );
    forwRate_sum    = zeros( nEvents, size(sampleStrike,2), length(sampleEpiPos) ); 
    forwRate_logl   = zeros( nEvents, size(sampleStrike,2), length(sampleEpiPos) );
    nAftershocks    = zeros( nEvents, size(sampleStrike,2), length(sampleEpiPos) );
    
    %% Computations via nested loops over sample and events
    % Evaluate model functions according to methodological settings
    ModelFuncs = set_modelFunctions( Inputs, ModelFuncs, sqrt( Inputs.ModelSettings.iniParamETAS ) );  
    
    % Apply time window 
    isInTime = t > t(idxAniso)' & t-t(idxAniso)' <= days(hours(Inputs.SpaceSettings.tStrikeEstim_hrs));                                                              
    
    % Loop over sample sets (nodal plane sample or epicenter position sample)
    for iStrike = 1:size(sampleStrike, 2)
        
        for iEpiPos = 1:length(sampleEpiPos)
                                                                                      
            % Compute start and end points of rupture segments
            [ rupStart, rupEnd ] = compute_rupCoords( x( idxAniso ), ...
                                                      y( idxAniso ), ...
                                                      rupExtent( idxAniso ), ...
                                                      sampleStrike( :, iStrike ), ...
                                                      sampleEpiPos( iEpiPos ) );        

            % Loop over all anisotropic events
            iAnisoEv = 0;
            for iEvent = idxAniso'

                iAnisoEv = iAnisoEv + 1;
                spatWindow = (Inputs.SpaceSettings.restrFactor + 0.5) * rupExtent(iEvent);
                        
                [ forwRate_sum(iAnisoEv, iStrike, iEpiPos), ...
                  forwRate_logl(iAnisoEv, iStrike, iEpiPos), ...
                  nAftershocks_matrix(iAnisoEv, iStrike, iEpiPos)] = estimate_forwardRateContrib( isInTime, ...
                                                                                                  iAnisoEv, ...
                                                                                                  x, y, mag, rupExtent, ...
                                                                                                  iEvent, ...
                                                                                                  spaceUnit, ...
                                                                                                  spatWindow, ...
                                                                                                  rupStart, rupEnd, ...
                                                                                                  minDist, ...
                                                                                                  ModelFuncs, ...
                                                                                                  q );
                                                                                          
            end              
        end
    end  

    %% Evaluation
    % Find sample indices that maximize distributed ETAS rates per event
%     isDuplicate = Catalog.evWeight(idxAniso)<1 & Catalog.id(idxAniso)==Catalog.id(idxAniso+1);
    [strikes, epiPos, nAftershocks] = find_maxTrigDensity( forwRate_sum, ...
                                                             nAftershocks_matrix, ...
                                                             sampleStrike, ...
                                                             sampleEpiPos, ...
                                                             Catalog.isDupl(idxAniso) );
                                                         
    %% Optional: Table with analysis results
    magAniso    = mag(idxAniso);
    tableAll    = sortrows(table(magAniso, nAftershocks, strikes, epiPos),'magAniso');
    tableAniso  = tableAll(tableAll.strikes>-1,:);
                                         
    %% Test isotropic kernel as well
    forwRateIso_sum = zeros(nEvents, 1);
    iAnisoEv        = 0;
    rupStart        = [x(idxAniso), y(idxAniso)];
    rupEnd          = rupStart;
    
    for iEvent = idxAniso'
        iAnisoEv    = iAnisoEv + 1;
        spatWindow  = Inputs.SpaceSettings.restrFactor * rupExtent(iEvent);
        forwRateIso_sum(iAnisoEv) = estimate_forwardRateContrib( isInTime, ...
                                                                  iAnisoEv, ...
                                                                  x, y, mag, rupExtent, ...
                                                                  iEvent, ...
                                                                  spaceUnit, ...
                                                                  spatWindow, ...
                                                                  rupStart, rupEnd, ...
                                                                  minDist, ...
                                                                  ModelFuncs, ...
                                                                  q );
                                                              
        if forwRateIso_sum(iAnisoEv) > max(max(forwRate_sum(iAnisoEv,:,:))) || nAftershocks(iAnisoEv) < 2
            strikes(iAnisoEv)   = -1;
            epiPos(iAnisoEv)    = -1;
        end
    end
    
    %% Optional: Plotting
    makePlot = false;
    if makePlot
%         close all
        idxEvents2plot  = idxAniso([1,3]);
        latFactor       = cos(35.769 * pi / 180);
        isInTime_1h     = t > t(idxEvents2plot)' & t-t(idxEvents2plot)' <= days(hours(1)); 
        isInTime_30h    = t > t(idxEvents2plot)' & t-t(idxEvents2plot)' <= days(hours(30));
        
        %% Plot aftershocks around estimated rupture line segment
        % Ridgecrest M6.4 event
        iEvent = idxEvents2plot(1);
        [ rupStart, rupEnd ] = compute_rupCoords( Catalog.x(iEvent), ...
                                                  Catalog.y(iEvent), ...
                                                  rupExtent(iEvent), ...
                                                  strikes(1:2), ...
                                                  epiPos(1:2) );
        figure
        s1 = scatter(Catalog.x(isInTime_30h(:,1)), Catalog.y(isInTime_30h(:,1)), 10, [0.9290 0.6940 0.1250], 'filled');
        hold on
        s2 = scatter(Catalog.x(isInTime_1h(:,1)), Catalog.y(isInTime_1h(:,1)), 20, 'r', 'filled');
        l1 = plot([rupStart(1,1), rupEnd(1,1)], [rupStart(1,2), rupEnd(1,2)], 'k-', 'LineWidth', 2);
        plot([rupStart(2,1), rupEnd(2,1)], [rupStart(2,2), rupEnd(2,2)], 'k-', 'LineWidth', 2);
        p1 = plot( Catalog.x(iEvent), Catalog.y(iEvent), 'p', 'Color', 'k', 'MarkerSize', 20 );
        p1.MarkerFaceColor = 'y';
%         p1fake = plot( 999, 999, 'p', 'Color', 'k', 'MarkerSize', 50 );
%         p1fake.MarkerFaceColor = 'y';
        legend([s1, s2, l1, p1], {'Aftershocks within 30 hours', 'Aftershocks within 1 hour', ...
                                      'Estimated M6.4 faults', 'M6.4 epicenter'})
        axis([-0.07, 0.14, -0.23, 0.11])
        yticks(-0.1688:0.1:0.0312)
        yticklabels(35.6:0.1:35.8)
        ylabel('Longitude (°)')
        xticks([-0.0023, 0.0788])
        xticklabels([-117.6, -117.5])
        xlabel('Latitude (°)')
        make_titleLeftCorner('(c)')
        
        % Ridgecrest M7.1 event
        iEvent = idxAniso(3);
        [ rupStart, rupEnd ] = compute_rupCoords( Catalog.x(iEvent), ...
                                                  Catalog.y(iEvent), ...
                                                  rupExtent(iEvent), ...
                                                  strikes(3), ...
                                                  epiPos(3) );
        
        figure
        s1 = scatter(Catalog.x(isInTime_30h(:,2)), Catalog.y(isInTime_30h(:,2)), 10, [0.5843 0.8157 0.9882], 'filled');
        hold on
        s2 = scatter(Catalog.x(isInTime_1h(:,2)), Catalog.y(isInTime_1h(:,2)), 20, 'b', 'filled');
        l1 = plot([rupStart(1), rupEnd(1)], [rupStart(2), rupEnd(2)], 'k-', 'LineWidth', 2);
        p1 = plot( Catalog.x(iEvent), Catalog.y(iEvent), 'h', 'Color', 'k', 'MarkerSize', 20 );
        p1.MarkerFaceColor = 'y';
        legend([s1, s2, l1, p1], {'Aftershocks within 30 hours', 'Aftershocks within 1 hour', ...
                                  'Estimated M7.1 fault', 'M7.1 epicenter'})
        axis([-0.35, 0.35, -0.35, 0.52])
        yticks(-0.2688:0.1:0.4312)
        yticklabels(35.5:0.1:36.2)
        ylabel('Longitude (°)')
        xticks([-0.3804, -0.2993, -0.2181, -0.1370, -0.0558, 0.0253, 0.1064, 0.1876, 0.2687, 0.3498])
        xticklabels(35.3:0.1:36.2)
        xlabel('Latitude (°)')
        make_titleLeftCorner('(d)')
                              
    end
                   
end
%% Old code pieces
%         plot_fittedStrike( Catalog, ...
%                               idx2plot, ...
%                               strikes2plot, ...
%                               epiPos2plot, ...
%                               x, y, t, mag, wi, rupExtent, ...
%                               Catalog.lon, Catalog.lat, ...
%                               spaceUnit, ...
%                               spatRestr, ...
%                               minSpaceDist, ...
%                               g_inn, f_inn, p, q )
%                           
%         plot_densityVsAftershocks( isInTime, ...
%                                    idx2plot, ...
%                                    strikes2plot, ...
%                                    epiPos2plot, ...
%                                    Catalog.x, Catalog.y, Catalog.mag, Catalog.rupExtent, ...
%                                    paramETAS, ...
%                                    spaceUnit, ...
%                                    f_inn, spatialFactor )
