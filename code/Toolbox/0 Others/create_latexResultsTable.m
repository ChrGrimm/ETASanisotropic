function create_latexResultsTable( tModelComparison, tableFormat )

    %% Format 1: Columns = Models; Rows = Variables
    % mu
    printResults( round(tModelComparison.mu, 2), '$\\mu$', '$\\frac{1}{days}$', '.2f' )
    % A
    printResults( round(tModelComparison.A, 3), '$A$', '', '.3f' )
    % alpha
    printResults( round(tModelComparison.alpha, 2), '$\\alpha$', '$\\frac{1}{mag}$', '.2f' )
    % c
    printResults( round(tModelComparison.c, 3), '$c$', '$\\frac{1}{days}$', '.3f' )
    % p
    printResults( round(tModelComparison.p, 2), '$p$', '', '.2f' )
    % D
    printResults( round(tModelComparison.D_km, 3), '$D$', '$km^2$', '.3f' )
    % gamma
    printResults( round(tModelComparison.gamma, 2), '$\\gamma$', '$\\frac{1}{mag}$', '.2f' )
    % q
    printResults( round(tModelComparison.q, 2), '$q$', '', '.2f' )
    % Tb
    printResults( round(tModelComparison.Tb_sek, 2), '$T_b$', '$sec$', '.1f' )
    % b
    printResults( round(tModelComparison.b, 2), '$b$', '', '.2f' )
    % Log-Likelihood
    printResults( round(tModelComparison.loglik_fval), '$l(\\theta)$', '', 'd' )
    % Branching Ratio
    printResults( round(tModelComparison.branchRatio1, 2), '$BR$', '', '.2f' )

    % Function definition
    function printResults( data, strVariable, strUnit, format )
        format      = ['%', format, ' & '];
        fprintf([strVariable, ' & ', strUnit, ' & ', format, format, format, format, format, '& ', format, format, format, format, format(1:end-2), ' \\\\ \n'], ...
                data(1), data(2), data(3), data(4), data(5), data(6), data(7), data(8), data(9), data(10))
    end

    %% Format 2: Columns = Variables; Rows = Models

end