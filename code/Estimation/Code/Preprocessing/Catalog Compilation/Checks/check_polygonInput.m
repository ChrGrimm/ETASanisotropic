function polygon = check_polygonInput( polygon )

    % Polygon numeric
    if any(any(isnan(polygon))) || any(any(~isnumeric(polygon)))
        error('NaN or non-numeric in ''polygon''')
    end
    % Polygon must be in counter-clockwise order
    if ispolycw(polygon(:,1),polygon(:,2))
        polygon = polygon(end:-1:1,:);
        disp('Polygon was made counter-clockwise')
    end
    % Polygon must be closed
    if ~all(polygon(end,:)==polygon(1,:))
        polygon(end+1,:) = polygon(1,:);
        disp('Polygon was closed by adding end-point')
    end
    % Polygon must have positive area
    if polyarea(polygon(:,1),polygon(:,2)) <= 0
        error('Polygon has no/ negative area')
    end
    % Polygon must have column dimension 2
    if ~( size(polygon,2)==2 )
        error('Polygon must have column dimension 2')
    end

end