function [ isInShadow, ...
           idxMainshock ] = check_ifMainshockBefore( Catalog, ...
                                                     iEvent, ...
                                                     idxInTime, ...
                                                     xCoord, ...
                                                     yCoord, ...
                                                     spaceUnit, ...
                                                     spatialWindow, ...
                                                     mag, ...
                                                     strMethod )

    if strcmp( strMethod, 'window' )
        
        %% Check if event is in the shadow of stronger "mainshock"
        % Check spatial shadow
        isInSpace   = getDistance( xCoord(idxInTime), yCoord(idxInTime), ...
                                   xCoord(iEvent),    yCoord(iEvent), spaceUnit ) ...
                          <= spatialWindow(idxInTime);
        % Check magnitude shadow
        isInMagn     = mag(idxInTime(isInSpace)) > mag(iEvent);
        isInShadow   = any( isInMagn );
        
        if isInShadow
            idxTimeAndSpace = idxInTime(isInSpace);
            idxMainshock    = idxTimeAndSpace(isInMagn);
        else
            idxMainshock    = [];
        end
        
        
    elseif strcmp( strMethod, 'cluster' )
        
        
        
    end

end