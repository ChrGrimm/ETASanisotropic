function [ r, ...
           r0, ...
           r_squ, ...
           r0_squ, ...
           segmFactor, ...
           isStartP, ...
           isEndP ] = compute_anisotropicGrid( polyStart, ...
                                               polyEnd, ...
                                               polyEdgeCategory, ...
                                               rupFlag, ...
                                               rupLength, ...
                                               rupXleft, ...
                                               rupXright, ...
                                               rupY, ...
                                               iEvent, ...
                                               spaceUnit, ...
                                               makeChecks, ...
                                               makePlots ) 
    
    %% Polygon data is empty, i.e. no anisotropic polygon part
    % Return empty vectors
    if isempty(polyStart)
        
        r           = [];
        r0          = [];
        r_squ       = [];
        r0_squ      = [];
        segmFactor  = [];
        isStartP    = true(0,1);
        isEndP      = true(0,1);
      
    %% Polygon data is not empty, i.e. there is anisotropic polygon part
    else
        
        %% Precomputations and initializations
        % Transpose start and end points of rotated polygon
        polyStart   = polyStart';
        polyEnd     = polyEnd';

        % Number of polygon edges
        nEdge   = size(polyStart,1);

        % Compute number of divisions per polygon edge
        if strcmp(spaceUnit, 'degree')
            divPerSpaceUnit = 50;
        elseif strcmp(spaceUnit, 'km')
            divPerSpaceUnit = 0.5;
        end
        nDiv        = max(1,ceil( divPerSpaceUnit*abs(polyEnd(:,2)-polyStart(:,2)) ));

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

        %% Compute coordinates x and y and booleans isStartP and isEndP
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

        %% Compute distances and weights of polygon grid segments
        % Compute distances of polygon grid points to rupture line as vertical difference since rupture
        % line was rotated to be horizontal
        % If distance = 0, add numerical zero.
        r           = abs( y-rupY );
        r(r<10^-15) = 10^-15;
        r_squ       = r.^2;

        % Compute third evaluation point along polygon grid segments
        y0 = y(isStartP) + r(isStartP)./(r(isStartP)+r(isEndP)) .* (y(isEndP)-y(isStartP));

        % Compute distance of third evaluation points to rupture line
        r0_squ                  = abs( y0-rupY ).^2;
        r0_squ(r0_squ<10^(-30)) = 10^(-30);
        r0                      = sqrt(r0_squ);

        % Compute weight of polygon grid segments with respect to entire anisotropic area only (!!!)
        segmWeight = abs(x(isEndP)-x(isStartP)) / (2*rupLength);

    %     % Neglect polygon grid edges that have very small distance to rupture line
    %     zero = abs(x(isEndP)-x(isStartP)).*r0_squ >= 10^-10 & r_squ(isStartP)+r_squ(isEndP) > 10^-20;

        % Create sign vector indicating whether polygon grid segment contributes positively or
        % negatively to spatial integral
        sign            = repelem(polyEdgeCategory', nDiv, 1);

        % Compute segment factor (positively/negatively oriented weight)
        segmFactor    = sign .* segmWeight; % .* zero;

        %% Checks
        if makeChecks

            try
                check_anisotropicGrid( segmWeight, ...
                                       segmFactor, ...
                                       iEvent, ...
                                       rupFlag, ...
                                       rupXleft, ...
                                       rupXright, ...
                                       rupY, ...
                                       rupLength, ...
                                       polyEdgeCategory, ...
                                       nDiv, ...
                                       dxx )
            catch
                warning('Fails check')
            end
        end
        
        %% Plots
        if makePlots

            plot_anisotropicGrid( x, ...
                                  y, ...
                                  r, ...
                                  rupXleft, ...
                                  rupXright, ...
                                  rupY, ...
                                  polyStart, ...
                                  polyEnd, ...
                                  polyEdgeCategory )
             
        end
        
    end
        
end