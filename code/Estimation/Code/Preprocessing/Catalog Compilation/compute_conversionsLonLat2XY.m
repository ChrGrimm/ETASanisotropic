function ModelFuncs = compute_conversionsLonLat2XY( ModelFuncs, polygonTarget )

    [lonCenter, latCenter]      = centroid( polyshape(polygonTarget(:,1), polygonTarget(:,2)) );
    latFactor                   = cos(latCenter * pi / 180);
    ModelFuncs.fConvert_lon2x   = @(lon) (lon-lonCenter) * latFactor;
    ModelFuncs.fConvert_lat2y   = @(lat) lat - latCenter;
    ModelFuncs.fConvert_x2lon   = @(x) x/latFactor + lonCenter;
    ModelFuncs.fConvert_y2lat   = @(y) y + latCenter;
    
end