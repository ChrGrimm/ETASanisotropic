function [dist2line, maxDist, lengthLine] = get_dist2line( linesStart, ...
                                                            linesEnd, ...
                                                            points, ...
                                                            spaceUnit, ...
                                                            isCalcDetails )
%
% This function computes the distance (and additional infos) between a set of lines given by
% their respective start and end coordinates and a set of points.
%
% Call: [dist2line, maxDist, lengthLine] = get_distanceLine( linesStart, linesEnd, points, addInfo )
%
% INPUTS:
%
% - linesStart: numeric array (Nx2)     x-y coordinates of N line start points
% - linesEnd:   numeric array (Nx2)     x-y coordinates of N line end points
% - points:     numeric array (2xM)     x-y coordinates of M points
% - addInfo:    boolean                 states whether additional info should be computed
%
% OUTPUTS:
%
% - dist2line:  numeric array (NxM)     distances between given N lines and M points
% - maxDist:    numeric array (NxM)     maximal distance between given N lines and M points 
% - lengthLine: numeric vector (Nx1)    length of line

    %% Initialization and dimension checks
    % Initialize distance matrix
    dist2line = zeros( size(linesStart,1), size(points,2) );
    % Checks
    if ~all( size(linesStart) == size(linesEnd) )
        error('LineStart and LineEnd have different dimensions!')
    elseif ~( size(linesStart,2)==2 )
        error('LineStart has not column dimension 2!')
    elseif ~( size(points,1)==2 )
        error('Point has not row dimension 2!')
    elseif size(points,2)>1 && size(linesStart,1)>1
        error('Either lineStart/End or points may be list of points')
    end
    
    isMultiPoints = size(points,2)>1;

    %% Directional vector and projection onto line
    % directional vector
    m = linesEnd - linesStart;

    % projection parameter
    proj = ( m(:,1) .* (points(1,:)-linesStart(:,1)) + ...
             m(:,2) .* (points(2,:)-linesStart(:,2)) ) ./ ( m(:,1).^2 + m(:,2).^2 );

    %% compute distance of points to line ...
    isProj0orNaN        = proj <= 0 | isnan(proj);
    isProj01            = proj >  0 & proj < 1;
    isProj1             = proj >= 1;
    
    % ... if projection to line is before linesStart or is actually a point (proj=NaN)
    if sum(isProj0orNaN)>0
        if isMultiPoints
            dist2line(isProj0orNaN) = getDistance( points(1,isProj0orNaN),  points(2,isProj0orNaN), ...
                                                   linesStart(:,1),         linesStart(:,2), ...
                                                   spaceUnit );
        else
            dist2line(isProj0orNaN) = getDistance( points(1,:),                 points(2,:), ...
                                                   linesStart(isProj0orNaN,1),  linesStart(isProj0orNaN,2), ...
                                                   spaceUnit );
        end
    end
        
    % ... if projection to line is after linesEnd
    if sum(isProj1)>0
        if isMultiPoints
            dist2line(isProj1) = getDistance( points(1,isProj1),    points(2,isProj1), ...
                                              linesEnd(:,1),        linesEnd(:,2), ...
                                              spaceUnit );
        else
            dist2line(isProj1) = getDistance( points(1,:),          points(2,:), ...
                                              linesEnd(isProj1,1),  linesEnd(isProj1,2), ...
                                              spaceUnit );
        end
    end
    
    % ... if projection to line is between linesStart and linesEnd
    if sum(isProj01)>0
        if isMultiPoints
            linePointLat        = linesStart(:,1) + proj(isProj01) .* m(:,1);
            linePointLon        = linesStart(:,2) + proj(isProj01) .* m(:,2);
            dist2line(isProj01) = getDistance( points(1,isProj01),  points(2,isProj01), ...
                                               linePointLat,        linePointLon, ...
                                               spaceUnit );
        else
            linePointLat        = linesStart(isProj01,1) + proj(isProj01) .* m(isProj01,1);
            linePointLon        = linesStart(isProj01,2) + proj(isProj01) .* m(isProj01,2);
            dist2line(isProj01) = getDistance( points(1,:),     points(2,:), ...
                                               linePointLat,     linePointLon, ...
                                               spaceUnit );
        end
    end
    
    if isCalcDetails
        % compute maximal distance to line
        maxDist = max( [ getDistance( points(1,:), points(2,:), linesStart(:,1), linesStart(:,2), spaceUnit ), ...
                         getDistance( points(1,:), points(2,:),  linesEnd(:,1), linesEnd(:,2), spaceUnit ) ], ...
                       [], 2 );

        % compute length of the line
        lengthLine = getDistance( linesStart(:,1), linesStart(:,2), linesEnd(:,1), linesEnd(:,2), spaceUnit );
    else
        maxDist     = 0;
        lengthLine  = 0;
    end
        
end