function [ distance, cdfValue ] = compute_spatialCDF( Mw, ...
                                                        Method, ...
                                                        f_inn, ...
                                                        q, ...
                                                        strDistance )

    %% Evaluate spatial cdf for given magnitude
    % Estimate rupture dimensions
    [rupL, ~] = estimate_rupSize( Mw, Tectonics.type, Tectonics.faultingStyle, Method.spaceUnit );
    % Define spatial extent
    k               = Method.factorSpatialRestr;
    c               = Method.minSpatialRestr;
    spatialExtent   = max(k * rupL, c);

    % Compute spatial cdf
    distance        = 0:rupL/1000:spatialExtent;
    if contains( Method.spatialDesign, 'normalized')
        normalizFactor = (1-f_inn(spatialExtent, spatialExtent.^2, Mw, rupL, 1).^(1-q));
    else
        normalizFactor = 1;
    end
    cdfValue        = (1-f_inn(distance, distance.^2, Mw, rupL, 1).^(1-q)) / normalizFactor;
    
    %% If opted for: Distance scaled by rupture length
    if strcmp( strDistance, 'scaled' )
        distance = distance/rupL;
    end
        
end