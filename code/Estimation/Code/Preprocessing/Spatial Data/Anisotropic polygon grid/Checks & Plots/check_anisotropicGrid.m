function check_anisotropicGrid( segmWeight, ...
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
    
    % Error tolerance 
    eps = 0.25 * 10^(-3);
    
    % Save rupture coordinates
    rupCoord            = [ rupXleft, rupY; rupXright, rupY ];                        
                            
    % Check if sum of weights and factors is between 0 and 1
    if sum(segmWeight) < 0 || sum(segmFactor) < -eps
        error(['Event ', num2str(iEvent), ': Sum of segmWeight or -segmFactor is smaller than 0!'])
    elseif sum(segmFactor) > 1+eps % || sum(segmWeight) > 1+eps 
        error(['Event ', num2str(iEvent), ': Sum of segmWeight or -segmFactor larger than 1!'])
    end

%     % Determine whether rupture start and end point are inside of polygon
%     isRuptureInPolygon  = inpolygon( rupCoord(:,1), rupCoord(:,2), polygon(1,:)', polygon(2,:)' );

    switch sum(rupFlag)
        % Case: Rupture is completely outside of polygon
        % NOTE: This check fails in the unlikely case that rupture end points are outside, but
        % inner parts of rupture line are inside of polygon
        case 0
            % Sum of weights must meet expected value, sum of factors must be 0
            if abs( sum(segmWeight) - sum(nDiv .* abs(dxx)) / (2 * rupLength) ) > eps
                error(['Event ', num2str(iEvent), ': Sum of segmWeight is unexpected value'])
            elseif abs(sum(segmFactor)) > eps
                error(['Event ', num2str(iEvent), ': Sum of segmFactor must be 0!'])
            end
        % Case: One end of rupture is inside of polygon, the other is outside of polygon
        case 1
            % Sum of weights and factors must meet expected values
            if abs( sum(segmWeight) - sum(nDiv .* abs(dxx)) / (2*rupLength) ) > eps
                error(['Event ', num2str(iEvent), ': Sum of segmWeight is unexpected value'])
            elseif abs( sum(segmFactor) - sum(nDiv .* abs(dxx) .* polyEdgeCategory') / (2 * rupLength) ) > eps
                error(['Event ', num2str(iEvent), ': Sum of segmFactor is unexpected value'])
            end
        % Case: Both end points of rupture are inside of polygon
        % NOTE: This check fails in non-rectangular polygon shapes
        case 2
%             if abs(sum(segmWeight) - 1) > eps
%                 error(['Event ', num2str(iEvent), ': Sum of segmWeight must be 1!'])
%             else
            if abs(sum(segmFactor) - 1) > eps
                error(['Event ', num2str(iEvent), ': Sum of segmFactor must be 1!'])
            end
    end
                            
end