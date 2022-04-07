function [ isInExtent, ...
           dist, ...
           dist2, ...
           rupStart, ...
           rupEnd ] = precompute_crossEventDistances( Catalog, spaceUnit, minXEventDistance )
       
    %% Initializations of new catalog columns 
    % Boolean: Marks which events are within spatio-temporal extent
    isInExtent  = cell(size(Catalog,1),1);
    % Number of events that are within spatio-temporal extent
    nContacts   = zeros(size(Catalog,1),1); 
    % Column of (squared) distances to all events within spatio-temporal extent
    dist        = cell(size(Catalog,1),1);
    dist2       = cell(size(Catalog,1),1);   
    % Start and end point of rupture line (for events modeled anisotropically)
    rupStart    = zeros(size(Catalog,1),2);
    rupEnd      = zeros(size(Catalog,1),2);
    
    %% Extract columns from table Catalog
    x           = Catalog.x;
    y           = Catalog.y;
    t           = Catalog.t;
    rupExtent   = Catalog.rupExtent;
    strike      = Catalog.strike;
    epiPos      = Catalog.epiPos;
    flag        = Catalog.flag;
    spatRestr   = Catalog.spatRestr;
    tempRestr   = Catalog.tempRestr;
    
    %% Computations
    % Loop over all events
    for iEvent = 1:size(Catalog,1)

        %% Compute distances to all past events
        % Compute start and end points of rupture lines
        % For events that are modeled isotropically, rupStart=rupEnd
        [ rupStart(iEvent,:), ...
          rupEnd(iEvent,:) ] = compute_rupCoords( x(iEvent), ...
                                                    y(iEvent), ...
                                                    rupExtent(iEvent), ...
                                                    strike(iEvent), ...
                                                    epiPos(iEvent) ); 
                                                
        % Compute distances between past and current event
        % For events that are modeled isotropically, the point-to-point distance is computed
        distSpace = get_dist2line( rupStart(1:iEvent-1,:), ...
                                   rupEnd(1:iEvent-1,:), ...
                                   [ x(iEvent), y(iEvent) ]', ...
                                   spaceUnit, ...
                                   0 );
           
        %% Only for target events: Store distances and evaluate events within spatio-temporal extent
        if flag(iEvent) > 0
            
            % Initialize boolean marker for event iEvent
            isInExtent{iEvent} = false(size(x));
            % Time difference in days to past events
            distTime = t(iEvent) - t(1:iEvent-1);
            % Evaluate boolean for past events (future events remain false in any case)
            isInExtent{iEvent}(1:iEvent-1) = distSpace <= spatRestr(1:iEvent-1) ... 
                                            & distTime >  0 ...
                                            & distTime <= tempRestr(1:iEvent-1);
            % Count past events that are within spatio-temporal extent
            nContacts       = nContacts + isInExtent{iEvent};
            % Store (squared) distances for events within spatio-temporal extent only
            dist{iEvent}    = max( distSpace(isInExtent{iEvent}(1:iEvent-1)), minXEventDistance );
%             dist{iEvent}    = max( minDist, distSpace( isInExtent{iEvent}(1:iEvent-1) ) );
            dist2{iEvent}   = dist{iEvent}.^2;
            
        end

    end
       
end