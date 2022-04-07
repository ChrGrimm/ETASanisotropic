function RotationMatrix = compute_rotationMatrix( strike, strDirection )

    if strcmp(strDirection, 'horizontal')
        alpha = (strike-90)/360 * 2*pi;
    elseif strcmp(strDirection, 'vertical')
        alpha = strike/360 * 2*pi;
    else
        error('Unknown rotation direction')
    end
        
    RotationMatrix = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
    
end