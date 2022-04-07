function [ rupStart, rupEnd ] = compute_rupCoords( epiX, epiY, rupL, strike, pos )

    if rupL == 0
        rupStart    = [epiX, epiY];
        rupEnd      = [epiX, epiY];
    else
        % Convert strike to mathematical rotation angle
        alpha       = 90-strike;
        % Compute rupture start and end point
        rupStart    = [epiX, epiY] - pos .* rupL .* [ cosd(alpha), sind(alpha) ];
        rupEnd      = [epiX, epiY] + (1-pos) .* rupL .* [ cosd(alpha), sind(alpha) ];
    end
    
end