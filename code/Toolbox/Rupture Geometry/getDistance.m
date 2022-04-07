function dist = getDistance( lon1, lat1, lon2, lat2, method )
%Returns geographical distance in km. lat1 and lon1 can be column vectors.
% The mean earth radius in km

    if strcmp(method, 'km')
        R = 6371.3;

        lat1 = deg2rad(lat1); lon1 = deg2rad(lon1);
        lat2 = deg2rad(lat2); lon2 = deg2rad(lon2);

        a = (sin((lat2-lat1)/2)).^2 + ((sin((lon2-lon1)/2)).^2).*cos(lat1).*cos(lat2);

        dist = 2*R*atan2(sqrt(a), sqrt(1-a));

    elseif strcmp(method, 'degree')

        dist = sqrt((lon1-lon2).^2 + (lat1-lat2).^2);

    else

        error('Invalid choice of space unit!')

    end

end
