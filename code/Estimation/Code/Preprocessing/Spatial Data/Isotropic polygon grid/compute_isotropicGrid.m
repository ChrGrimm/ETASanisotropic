function [ r, ...
           r0, ...
           r_squ, ...
           r0_squ, ...
           segmFactor, ...
           isStartP, ...
           isEndP ] = compute_isotropicGrid( polyStart, ...
                                             polyEnd, ...
                                             iEvent, ...
                                             refX, ...
                                             refY, ...
                                             refFlag, ...
                                             spaceUnit, ...
                                             version, ...
                                             makeChecks, ...
                                             makePlots ) 

    %% Return empty vectors, if polygon data is empty
    % i.e. no isotropic/isotropic-left/isotropic-right polygon part
    if isempty(polyStart)
        
        r           = [];
        r0          = [];
        r_squ       = [];
        r0_squ      = [];
        segmFactor  = [];
        isStartP    = true(0,1);
        isEndP      = true(0,1);
      
    %% Compute polygon-distance data if polygon data is not empty
    % i.e. there is isotropic/isotropic-left/isotropic-right polygon part
    else

        %% Prepare polygon grid
        % Number of polygon edges
        nEdge   = size(polyStart,1);
        
        % Compute length of polygone edges, and minimum/maximum distance of location to edge
        [dist2line, maxDist, lengthLine] = get_dist2line(polyStart, polyEnd, [refX; refY], spaceUnit, 1); 

        % Compute number of divisions per polygon edge
        if strcmp(spaceUnit, 'degree')
            divPerSpaceUnit = 50;
        elseif strcmp(spaceUnit, 'km')
            divPerSpaceUnit = 0.5;
        end
        nDiv = max( 1, ceil( divPerSpaceUnit*lengthLine .* (1 - dist2line ./ maxDist) ) );

        % Determine x and y increment by division of each polygon edge
        dxx = ( polyEnd(:,1)-polyStart(:,1) ) ./ nDiv;
        dyy = ( polyEnd(:,2)-polyStart(:,2) ) ./ nDiv;

        % Initialize polygon grid coordinate vectors x and y (total number of divisions + one starting
        % grid point per edge)
        x = zeros( sum(nDiv)+nEdge, 1 ); 
        y = x;

        % Initialize boolean vector indicating whether x-y pair is start and end point of a polygon grid
        % segment
        isStartP    = true(size(x));
        isEndP      = isStartP;

        %% Evaluate polygon grid coordinates x and y and booleans isStartP and isEndP
        for iEdge = 1:nEdge

            % indices of start and end point of current polygon edge
            idxStartP               = sum(nDiv(1:iEdge-1)) + iEdge;
            idxEndP                 = sum(nDiv(1:iEdge))   + iEdge;

            % Only end point of current polygon edge is not a grid segment start point
            isStartP( idxEndP )     = false;
            % Only start point of current polygon edge is not a grid segment end point
            isEndP( idxStartP )     = false;

            % Fill according ranges of x and y with polygon grid point coordinates
            x( idxStartP:idxEndP )  = polyStart(iEdge,1) + dxx(iEdge) * (0:nDiv(iEdge));
            y( idxStartP:idxEndP )  = polyStart(iEdge,2) + dyy(iEdge) * (0:nDiv(iEdge));

        end

        %% Compute distances of polygon grid coordinates to reference point
        % Distances between polygon grid coordinates and reference point
        r               = max( 10^-15, getDistance( x, y, refX, refY, spaceUnit ));
        r_squ           = r.^2;   
        % Distances between neighbouring grid coordinates forming a triangle with reference point
        r12_squ         = getDistance( x(isStartP), y(isStartP), x(isEndP), y(isEndP), spaceUnit ).^2;
               
        % Compute third polygon grid evaluation points (between neighboring grid coordinates)
        x0 = x(isStartP) + r(isStartP)./(r(isStartP)+r(isEndP)) .* (x(isEndP)-x(isStartP));
        y0 = y(isStartP) + r(isStartP)./(r(isStartP)+r(isEndP)) .* (y(isEndP)-y(isStartP));

        % Compute distance of third polygon grid evaluation points to reference point
        r0              = max(10^-15, getDistance( x0, y0, refX, refY, spaceUnit ));
        r0_squ          = r0.^2;
        
        %% Compute weights of segments by "Kosinussatz"
        % Compute radiant between neighbouring polygon grid coordinates and reference point
        % Note: Radiant computation only works properly with space unit 'degree', therefore 
        % distances between polygon grid coordinates and reference point need to be recomputed if 
        % space unit 'km' was used before. 
        if strcmp(spaceUnit, 'degree')
            radiant         = ( r_squ(isStartP) + r_squ(isEndP) - r12_squ ) ./ ( 2 * r(isStartP) .* r(isEndP) );
        elseif strcmp(spaceUnit, 'km')
            % Distances between polygon grid coordinates and reference point
            r_degree        = max(10^-15, getDistance( x, y, refX, refY, 'degree' ));
            r_degree_squ    = r_degree.^2; 
            % Distances between neighbouring grid coordinates forming a triangle with reference point
            r12_degree_squ  = getDistance( x(isStartP), y(isStartP), x(isEndP), y(isEndP), 'degree' ).^2;            
            radiant         = ( r_degree_squ(isStartP) + r_degree_squ(isEndP) - r12_degree_squ )...
                              ./  ( 2 * r_degree(isStartP) .* r_degree(isEndP) );   
        end
        % Correct radiant if necessary
        if any(abs(radiant)>1)
            radiant(abs(radiant)>1) = 1-10^-10;
            disp(['Set radiant to 1-10^-10 for iEvent ', num2str(iEvent), ' because abs(radiant)>1'])
        end        
        % Apply "Kosinussatz" to compute weight of segments spanned by triangle between neighbouring 
        % polygon grid coordinates and reference point (scaled by 360° circle)
        segmWeight              = acos(radiant) /(2*pi);  

        %% Compute segment factors by accounting for triangle orientation
        % Compute sign vector indicating whether polygon grid segment contributes positively or
        % negatively to spatial integral
        triangleDeterminant = (x(isStartP).*y(isEndP) + refX*y(isStartP) + refY*x(isEndP)) ...
                              - (x(isEndP).*y(isStartP) + refX*y(isEndP) + refY*x(isStartP));
        sign                = (-1) * (triangleDeterminant < 0) + 1 * (triangleDeterminant >= 0);

        % Neglect polygon grid edges that have very small distance to location or very small absolute 
        % value of trianglular determinant (= area of triangle)
        zero        = abs( triangleDeterminant ) >= 10^-12  &  r_squ(isStartP) + r_squ(isEndP) > 10^-20;

        % Compute segment factor (positively/negatively oriented weight)
        segmFactor  = sign .* segmWeight .* zero;
        
        %% Update polygon category for negative orientations
        polyEdgeCategory = sign(cumsum(nDiv))';
        
        %% Optional: Checks and plots 
        if makeChecks
            check_isotropicGrid( iEvent, refFlag, segmWeight, segmFactor, version )            
        end


        if makePlots
            
            plot_isotropicGrid( x, ...
                                y, ...
                                r, ...
                                refX, ...
                                refY, ...
                                polyStart, ...
                                polyEnd, ...
                                segmWeight, ...
                                segmFactor, ...
                                polyEdgeCategory, ...
                                version )
            
        end
        
    end

end