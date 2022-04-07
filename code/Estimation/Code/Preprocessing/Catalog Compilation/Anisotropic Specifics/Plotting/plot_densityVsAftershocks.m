function plot_densityVsAftershocks( isInTime, ...
                                   idx, ...
                                   strikes, ...
                                   epiPos, ...
                                   x, y, mag, rupL, ...
                                   paramETAS, ...
                                   spaceUnit, ...
                                   f_inn, spatialFactor )
                 
    isInTime = isInTime(:,1) | isInTime(:,2);
    
    %% Anisotropic Model
    % Compute Density at space grid
%     load('Pol_Ridgecrest_75km.mat')
    xGrid = min(x):0.005:max(x);
    yGrid = min(y):0.005:max(y);
    [xGrid,yGrid] = meshgrid(xGrid,yGrid);
    xGridVector = xGrid(:);
    yGridVector = yGrid(:);
    [ rupStart, rupEnd ] = compute_rupCoords( x(idx(1)), y(idx(1)), rupL(idx(1)), strikes{1}(1), epiPos{1}(1) );
    dist1 = reshape( get_dist2line( rupStart, rupEnd, [xGridVector'; yGridVector'], spaceUnit, 0 ), size(xGrid) );
    [ rupStart, rupEnd ] = compute_rupCoords( x(idx(1)), y(idx(1)), rupL(idx(1)), strikes{1}(2), epiPos{1}(2) );
    dist2 = reshape( get_dist2line( rupStart, rupEnd, [xGridVector'; yGridVector'], spaceUnit, 0 ), size(xGrid) );
    [ rupStart, rupEnd ] = compute_rupCoords( x(idx(2)), y(idx(2)), rupL(idx(2)), strikes{2}, epiPos{2} );
    dist3 = reshape( get_dist2line( rupStart, rupEnd, [xGridVector'; yGridVector'], spaceUnit, 0 ), size(xGrid) );
    
    A = paramETAS(2);
    alpha = paramETAS(3);
    q = paramETAS(8);
    productivity = A * exp(alpha*(mag(idx)-2.1));
    spatialFactor = [spatialFactor(2.5*rupL(idx(1)), mag(idx(1)), rupL(idx(1)), 1), ...
                     spatialFactor(2.5*rupL(idx(2)), mag(idx(1)), rupL(idx(2)), 1)];
    rate1 = productivity(1)/2 * spatialFactor(1) * f_inn(dist1, dist1.^2, mag(idx(1)), rupL(idx(1)), 1).^(-q);
    rate2 = productivity(1)/2 * spatialFactor(1) * f_inn(dist2, dist2.^2, mag(idx(1)), rupL(idx(1)), 1).^(-q);
    rate3 = productivity(2) * spatialFactor(2) * f_inn(dist3, dist3.^2, mag(idx(1)), rupL(idx(1)), 1).^(-q);
    
    % Plot Anisotropic
    figure
    scatter(xGridVector, yGridVector, 5, log10(rate1(:)+rate2(:)+rate3(:)), 'filled');
    colorbar
    hold on
    scatter( x(isInTime), y(isInTime), 1, 'b', 'filled' )
    xlabel('Longitude (centralized)', 'FontSIze', 12)
    ylabel('Latitude (centralized)', 'FontSIze', 12)
    xlim([min(xGridVector), max(xGridVector)])
    ylim([min(yGridVector), max(yGridVector)])
    
    %% Isotropic Model
    dist1 = reshape( getDistance( x(idx(1)), y(idx(1)), xGridVector, yGridVector, spaceUnit ), size(xGrid) );
    dist2 = reshape( getDistance( x(idx(2)), y(idx(2)), xGridVector, yGridVector, spaceUnit ), size(xGrid) );
    
    rate1 = productivity(1) * spatialFactor(1) * f_inn(dist1, dist1.^2, mag(idx(1)), 0, 1).^(-q);
    rate2 = productivity(2) * spatialFactor(2) * f_inn(dist2, dist2.^2, mag(idx(2)), 0, 1).^(-q);
    
    % Plot Anisotropic
    figure
    scatter(xGridVector, yGridVector, 5, log10(rate1(:)+rate2(:)), 'filled');
    colorbar
    hold on
%     isInTime = isInTime(:,1) | isInTime(:,2);
    scatter( x(isInTime), y(isInTime), 1, 'b', 'filled' )
    xlabel('Longitude (centralized)', 'FontSIze', 12)
    ylabel('Latitude (centralized)', 'FontSIze', 12)
    xlim([min(xGridVector), max(xGridVector)])
    ylim([min(yGridVector), max(yGridVector)])
    
end