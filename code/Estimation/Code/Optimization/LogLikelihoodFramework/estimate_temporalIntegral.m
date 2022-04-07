function [ G, dG_c, dG_p ] = estimate_temporalIntegral( inclGradients, ...
                                                        Catalog, ...
                                                        p, ...
                                                        fTempK_inn, ...
                                                        idx, ...
                                                        timeWindow )

    %% ETAS run with spatial parameter estimation
    % Evaluate inner temporal function ...
    % ... at lower integral boundary
    lowerBound  = min( max(0, timeWindow(:,1)'-Catalog.t(idx)), Catalog.tempRestr(idx) );
    g_inn_start = fTempK_inn( lowerBound );
    % ... at upper integral boundary
    upperBound  = max( 0, min( timeWindow(:,2)'-Catalog.t(idx), Catalog.tempRestr(idx) ) );
    g_inn_end   = fTempK_inn( upperBound );

    % Compute inner temporal results to the 1-p
    g_inn_start_1minusP = g_inn_start.^(1-p);
    g_inn_end_1minusP   = g_inn_end.^(1-p);

    % Compute temporal integral
    G          = 1/(1-p) * (g_inn_end_1minusP - g_inn_start_1minusP); % ./ tempScalingFactor;
    G          = sum(G,2);

    if inclGradients
        % Compute gradiants by c and p
        dG_c = g_inn_end.^(-p) - g_inn_start.^(-p);
        dG_p = 1/(1-p)^2 * (g_inn_end_1minusP - g_inn_start_1minusP) ...
               + 1/(1-p) * (-log(g_inn_end) .* g_inn_end_1minusP + log(g_inn_start) .* g_inn_start_1minusP);
        dG_c = sum(dG_c,2);
        dG_p = sum(dG_p,2);
    else
        dG_c = 0;
        dG_p = 0;
    end
        
end
