function [lon_shift, lat_shift] = shift_lonLatAlongStrike( lat, dist_km, strike )

    % mean latitude factor
    mLF = cosd(lat);   

    % coordinates shifted by dist_km along strike
    lon_shift = dist_km .* sind(strike) / mLF / 111;
    lat_shift = dist_km .* cosd(strike) / 111;
    
end
