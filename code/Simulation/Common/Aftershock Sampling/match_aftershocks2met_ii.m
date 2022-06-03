function Output = match_sampledLocation2met( x, y, MET, MET_locations )

    %%
    nAftershocks = length(x);
    outputMatrix = zeros(nAftershocks, 7);
    Output = array2table(outputMatrix, 'VariableNames',{'id', 'x', 'y', 'depth', 'mag', 'strike', 'epiPos'});
    
    maxMatrixDimension  = 2*10^8;
    nMET_locations      = size(MET_locations,1);
    nAftershocksAtOnce  = floor( maxMatrixDimension/nMET_locations );
    nIterations         = ceil( nAftershocks/nAftershocksAtOnce );
    
    locIDs              = zeros(nAftershocks, 1);
    
    for iIter = 1:nIterations
        idxAftershocks  = (iIter-1)*nAftershocksAtOnce + (1:nAftershocksAtOnce)';
        idxAftershocks(idxAftershocks>nAftershocks) = [];
        distances_km    = getDistance( x(idxAftershocks), y(idxAftershocks), ...
                                       MET_locations.x', MET_locations.y', 'km' );
        [~, locIDs(idxAftershocks)] = min(distances_km, [], 2); 
    end 
    
    uniqueLocIDs = unique(locIDs);
    
    for iLoc = uniqueLocIDs'
        isAftershInLoc  = locIDs == iLoc;
        nAftershInLoc   = sum(isAftershInLoc);
        idxEventsInLoc  = find(MET.locID == iLoc);
        if length(idxEventsInLoc)==1
            idxRand = idxEventsInLoc*ones(nAftershInLoc,1);
        else
            idxRand = randsample( idxEventsInLoc, nAftershInLoc, true, MET.eventProb(idxEventsInLoc) );
        end
        try
        Output(isAftershInLoc,:) = MET(idxRand,[1,5,6,9:12]);
        catch
            test=-1;
        end
        if mod(iLoc,250)==0
            disp(num2str(iLoc))
        end
    end

end