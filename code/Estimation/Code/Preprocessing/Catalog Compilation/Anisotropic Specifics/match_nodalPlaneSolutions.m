function tempCat = match_nodalPlaneSolutions( tempCat, polygon, strNodalPlaneSource )

    %% Load nodal plane source data
    if strcmp(strNodalPlaneSource, 'gcmt')
        load('nodalPlanes_gcmt_1976_2017.mat','NodalPlanes')
    elseif strcmp(strNodalPlaneSource, 'iscgem')
        load('nodalPlanes_iscGem_1905_2016.mat','NodalPlanes');
    elseif strcmp(strNodalPlaneSource, 'self') || strcmp(strNodalPlaneSource, 'none')
        return
    else
        error('No valid input in ''match_nodalPlaneSolutions.m''')
    end

    inpol       = inpolygon(NodalPlanes.lon, NodalPlanes.lat, polygon(:,1), polygon(:,2));
    NodalPlanes = NodalPlanes( inpol, : );
    
    %% Extract event time data from catalog
    nEvents = size(tempCat,1);

    %% Initialize strike-dip-rake columns
    tempCat.str1 = zeros(nEvents,1);
    tempCat.dip1 = zeros(nEvents,1);
    tempCat.rak1 = zeros(nEvents,1);
    tempCat.str2 = zeros(nEvents,1);
    tempCat.dip2 = zeros(nEvents,1);
    tempCat.rak2 = zeros(nEvents,1);

    %% Search matches
    vectorNoMatch       = [];
    counter             = zeros(3,1);
    
    for iNP = 1:size(NodalPlanes,1)
        % Time filter
        isTime = abs( minutes( tempCat.date - NodalPlanes.date(iNP) ) ) <= 2;
        % Location filter
        isLoc  = getDistance(tempCat.lon, tempCat.lat, NodalPlanes.lon(iNP), NodalPlanes.lat(iNP), 'km') ...
                    <= 2.5 * estimate_rupSize( NodalPlanes.mag(iNP), 'continental', 'all', 'km' ); % TO EDIT
        % Magnitude filter
        isMag  = abs( tempCat.mag - NodalPlanes.mag(iNP) ) <= 0.5;
        
        isMatch = isTime & isLoc & isMag;

        % If more than one event fulfills criteria: Select closes in space
        if sum( isMatch ) == 0

            vectorNoMatch = [vectorNoMatch; iNP];
            counter(1)  = counter(1)+1;

        else
            
            if sum(isMatch)>1
                counter(2) = counter(2)-1;
                counter(3) = counter(3)+1;
                distances = getDistance( tempCat.lon, tempCat.lat, ...
                                            NodalPlanes.longitude(iNP), NodalPlanes.latitude(iNP), 'degree' );
                [~,isMatch] = min(distances);
            end
            counter(2) = counter(2)+1;

            tempCat.str1(isMatch)   = NodalPlanes.strike1( iNP );
            tempCat.dip1(isMatch)   = NodalPlanes.dip1( iNP );
            tempCat.rak1(isMatch)   = NodalPlanes.rake1( iNP );
            tempCat.str2(isMatch)   = NodalPlanes.strike2( iNP );
            tempCat.dip2(isMatch)   = NodalPlanes.dip2( iNP );
            tempCat.rak2(isMatch)   = NodalPlanes.rake2( iNP );
            
        end

    end
    
end
