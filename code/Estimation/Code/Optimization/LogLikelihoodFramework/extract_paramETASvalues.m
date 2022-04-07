function [mu, A, alpha, c, p, D, gamma, q, Tb, beta] = extract_paramETASvalues( paramETAS, strMode )

    % Convert parameters to actual scale
    if strcmp(strMode, 'sqrt')
        paramETAS = paramETAS.^2;
    elseif ~strcmp(strMode, 'default')
        error('Unknown mode!')
    end
    
    % Extract parameters from vector
    mu      = paramETAS(1);
    A       = paramETAS(2);
    alpha   = paramETAS(3);
    c       = paramETAS(4);
    p       = paramETAS(5);
    D       = paramETAS(6);
    gamma   = paramETAS(7);
    q       = paramETAS(8);
    
    if length(paramETAS)==10
        Tb      = paramETAS(9);
        beta    = paramETAS(10);
    else
       Tb       = NaN;
       beta     = NaN;
    end

end