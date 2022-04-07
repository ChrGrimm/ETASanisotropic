function [ bandwidth, ...
           gaussDensity ] = precompute_backgrSpatDistr( Catalog, ...
                                                        spaceUnit, ...
                                                        spaceModel )
% This function precomputes the bandwidth and the corresponding isotropic Gaussian kernel for each
% event which is later used to redistribute background event locations.
%
% Call: [bwd, dGauss] = precompute_backgrSpatDistr( Catalog, spaceUnit )
%
% INPUTS:
%
% - Catalog:    table       historic event set. Fields used here:
%                           - flag: code for the kind of event with respect to target settings
%                           - x: (Scaled) longitude coordinates
%                           - y: (Scaled) latitude coordinates
% - spaceUnit:  string      describes the chosen space format (e.g. degree, kilometers, etc.)
%
% OUTPUTS:
%
% - bwd:        numeric vector (Nx1)    bandwidth per event
% - dGauss:     

    if strcmp( spaceModel, 'none' )
        
        bandwidth       = zeros(sum(Catalog.flag>=0),1);
        gaussDensity    = 1;
        
    else
        
        %% Bandwidth settings
        % Bandwidth defined as k-th nearest neighbour
        kNearest = 5;    
        % minimum bandwidth
        if strcmp(spaceUnit, 'km')
            minBandwidth = 6371.3 * pi / 180 * 0.05;
        elseif strcmp(spaceUnit, 'degree')
            minBandwidth = 0.05;
        end

        %% Compute bandwidth
        % Booleans: marks target time events, excluding duplicates of the same event
        isTargetT   = Catalog.flag >= 0 & ~Catalog.isDupl;
        % Compute distances between target time events
        distances        = getDistance( Catalog.x(isTargetT), Catalog.y(isTargetT), ...
                                        Catalog.x(isTargetT)', Catalog.y(isTargetT)', spaceUnit);
        % Sort distances
        distances_sorted = sort(distances, 2);
        % Define bandwidth as k-th largest distances (without distance of an event to itself)
        bandwidth        = max( distances_sorted(:, 1 + kNearest), minBandwidth );

        %% Compute Gaussian density
        % Extract distances from target time events to target time-space events
        distances       = distances(:, Catalog.flag(isTargetT)>0);
        bandwidth_squ   = bandwidth.^2;
        % Compute cross-event contributions according to Gaussian kernel 
        gaussDensity    = 1./(2*pi*bandwidth_squ) .* exp((-distances.^2)./(2*bandwidth_squ));
           
    end
    
end

