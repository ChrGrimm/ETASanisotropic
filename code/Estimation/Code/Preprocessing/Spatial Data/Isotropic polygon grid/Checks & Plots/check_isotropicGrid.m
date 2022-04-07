function check_isotropicGrid( iEvent, locFlag, segmWeight, segmFactor, version )

    % Error tolerance
    eps = 0.5 * 10^(-3);

    % Check if sum of weights and factors is between 0 and 1
    if sum(segmWeight) < -eps || sum(segmFactor) < -eps
        error(['Event ', num2str(iEvent), ': Sum of segmWeight or -segmFactor (', version, ') is smaller than 0!'])
    elseif sum(segmFactor) > 1+eps % || sum(segmWeight) > 1+eps
        error(['Event ', num2str(iEvent), ': Sum of segmWeight or -segmFactor (', version, ') is larger than 1!'])
    end

    switch abs( locFlag )
        % Case: Point lies outside of polygon
        case {0, 0.1}
            % Sum of weights must meet expected value, sum of factors must be 0
            if abs( sum(segmFactor) ) > eps
                error(['Event ', num2str(iEvent), ': Sum of segmFactor (', version, ') must be 0!'])
            end
        % Case: Point lies on polygon edge
        % NOTE: This check fails in non-rectangular polygon shapes
        case 0.5
            % Sum of weights and factors must meet expected values
            if strcmp(version, 'isotropic left') || strcmp(version, 'isotropic right')
                if ~(sum(segmWeight) <= 0.5+eps) || ~(sum(segmFactor) <= 0.5) 
                    error(['Event ', num2str(iEvent), ': Sum of segmWeight or -segmFactor (', version, ') must be smaller/equal to 0.5!'])
                end
            else
                if ~(sum(segmWeight) <= 1+eps) || ~(sum(segmFactor) <= 0.5+eps)
                    error(['Event ', num2str(iEvent), ': Sum of segmWeight or -segmFactor must be smaller/equal to 0.5!'])
                end
            end
        % Case: Point lies inside or on the edge of the polygon
        % NOTE: This check fails in non-rectangular polygon shapes
        case 1
            % Sum of weights and factors must meet expected values
            if strcmp(version, 'isotropic left') || strcmp(version, 'isotropic right')
                if abs(sum(segmFactor) - 0.5) > eps % || abs(sum(segmWeight) - 0.5) > eps
                    error(['Event ', num2str(iEvent), ': Sum of segmWeight or -segmFactor (', version, ') must be 0.5!'])
                end
            else
                if abs(sum(segmFactor) - 1) > eps % || abs(sum(segmWeight) - 1) > eps
                    warning(['Event ', num2str(iEvent), ': Sum of segmWeight or -segmFactor (', version, ') must be 1!'])
                end
            end
    end
            
end