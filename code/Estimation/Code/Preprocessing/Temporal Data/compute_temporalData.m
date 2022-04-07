function [ idxRecentEvents, ...
           idxTimeGrid, ...
           timeGrid, ...
           timeSpan4recentEvents, ...
           isZeroSpatialGradient, ...
           idxRecentEvents_spatGrad, ...
           idxTimeGrid_spatGrad ] = compute_temporalData( tFixGrid, ...
                                                            thresh4recentEvents, ...
                                                            thresh4largeMag, ...
                                                            c, p, ...
                                                            tEvents, ...
                                                            magEvents, ...
                                                            typeKernel, ...
                                                            timeGrid_interEvent, ...
                                                            maxTimeRestr_days )
   
    %% Initializations
    % Initialize output vectors
    nPreallocate            = 10^8;
    nTimeIntervals          = length(tFixGrid)/2; 
    timeGrid                = zeros( nPreallocate, 1 );
    idxRecentEvents         = zeros( nPreallocate, 1 );
    idxTimeGrid             = zeros( nPreallocate/10, 1 );
    timeSpan4recentEvents   = zeros( nTimeIntervals, 1 );
    % Initialize newIdx_data
    lastIdx_timeGrid        = 0;
    lastIdx_data            = 0;
    
    %% Loop over all time intervals
    for i = 1:nTimeIntervals
        % Determine start and end point of current time interval
        tStart                      = tFixGrid( 2*i-1 );
        tEnd                        = tFixGrid( 2*i );
        tLength                     = tEnd-tStart;
        % mean_interEventTime; % 0.5; % % 100 or 160
        timeSpan4recentEvents(i)    = max(0, nthroot(thresh4recentEvents./(p*(p+1)*tLength), -p-2) -c); 
        
        % Find contributing events
        % - Events that occur at most timeSpan4recentEvents(i) before current tStart
        % - EVents with magnitude > 5
        idxContribEvents            = find( tEvents <= tStart & ...
                                            tEvents >= tStart - maxTimeRestr_days & ...
                                            (tEvents >= tStart-timeSpan4recentEvents(i) | magEvents > thresh4largeMag) );
        
        % Determine time grid within current time interval
        interEvTimes                = tStart + timeGrid_interEvent( timeGrid_interEvent <= tLength );
        newTimeGrid                 = [ tStart, interEvTimes, tEnd ];
        nNewTimePoints              = length( newTimeGrid );
        newIdx_timeGrid             = lastIdx_timeGrid + ( 1:nNewTimePoints );
        lastIdx_timeGrid            = newIdx_timeGrid(end);
        timeGrid( newIdx_timeGrid ) = newTimeGrid;
        
        % Update output vectors
        nContribEvents                  = length( idxContribEvents );
        newIdx_data                     = lastIdx_data + ( 1:nContribEvents*nNewTimePoints );
        lastIdx_data                    = lastIdx_data + nContribEvents * nNewTimePoints;
        idxRecentEvents( newIdx_data )  = repmat( idxContribEvents, nNewTimePoints, 1 );
        idxTimeGrid( newIdx_data )      = repelem( newIdx_timeGrid, nContribEvents*ones(nNewTimePoints,1) );

    end

    % Delete unused (excess) rows of vectors
    deleteFromIdx                       = lastIdx_data + 1;
    idxRecentEvents( deleteFromIdx:end )= [];
    idxTimeGrid( deleteFromIdx:end )    = [];
    deleteFromIdx                       = lastIdx_timeGrid + 1;
    timeGrid( deleteFromIdx:end )       = [];
    
    %% Post-process results    
    isZeroSpatialGradient               = strcmp( typeKernel(idxRecentEvents), 'full' );
    idxRecentEvents_spatGrad            = idxRecentEvents( ~isZeroSpatialGradient );
    idxTimeGrid_spatGrad                = idxTimeGrid( ~isZeroSpatialGradient );

%     isZero_time                 = ismember( unique(idxTimeGrid), idxTimeGrid(~isZero_events) );
                            
end