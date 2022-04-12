function TriggerRelations = evaluate_triggRelationsZBZ( Catalog )

    %% 
    TriggerRelations = initialize_tTriggerRelations( Catalog, 'ZBZ' );
    
    %% Zaliapin Ben Zion Parameters
    % According to section 3.7 of Zaliapin, Ben-Zion (2016)
    b                   = 1.0; % log(beta)
    d                   = 1.3; % 1.3 (2016), 1.6 (2013)
    q                   = 0.5;
    dist_thresh         = 10^-3;

    %% Find nearest neighbors
    for iEvent = 2:size(Catalog,1)
        %% Compute distance metric
        % Compute temporal distance
        tij = Catalog.t(iEvent) - Catalog.t(1:iEvent-1);
        % Compute spatial distance
        rij = max( 0.1, getDistance( Catalog.lon(iEvent), Catalog.lat(iEvent), ...
                                     Catalog.lon(1:iEvent-1), Catalog.lat(1:iEvent-1), 'km' ) );
        % Compute distance measure including magnitude term
        mi  = Catalog.mag(1:iEvent-1);
        nij = tij .* rij.^d .* 10.^(-b*mi);       
        
        %% Evaluate closest distance
        [ minDistance, idxNearest ]             = min(nij); 
        if minDistance <= dist_thresh
            TriggerRelations.parentID(iEvent)   = Catalog.id(idxNearest);
        end
        TriggerRelations.distance_ZBZ(iEvent)   = minDistance;
        mi_nearest                              = mi(idxNearest);
        TriggerRelations.tij_rescaled(iEvent)   = tij(idxNearest)*10.^(-q*b*mi_nearest);
        TriggerRelations.rij_rescaled(iEvent)   = rij(idxNearest)^d*10.^(-(1-q)*b*mi_nearest);
        
    end
    
    %% Plots
    % Results ZBZ
    figure
    subplot(1,2,1)
    scatter(log10(TriggerRelationsZBZ.tij_rescaled), log10(TriggerRelationsZBZ.rij_rescaled), 5, 'b', 'filled')
    subplot(1,2,2)
    autoplot_smoothedHeatMap([log10(TriggerRelationsZBZ.tij_rescaled), log10(TriggerRelationsZBZ.rij_rescaled)])

    
end