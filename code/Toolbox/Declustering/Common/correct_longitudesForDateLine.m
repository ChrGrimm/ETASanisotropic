function longitudes = correct_longitudesForDateLine( longitudes )

    %% Convert all longitudes to positive values in order to avoid dateline issues
    longitudes = longitudes .* (longitudes>=0) + (360+longitudes) .* (longitudes<0);
    
end