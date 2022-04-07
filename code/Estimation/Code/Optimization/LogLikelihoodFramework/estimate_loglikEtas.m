function [ fval_logL, ...
            grad_logL, ...
            summand2_logL, ...
            anisoPeaks] = estimate_loglikEtas( inclGrad, ...
                                                 Catalog, ...
                                                 Inputs, ...
                                                 ModelFuncs, ...
                                                 sqrtParamETAS, ...
                                                 aggregBackgrRate, ...
                                                 SpatData )
                                             
    if any(sqrtParamETAS<0)
        test=-1;
    end
                                                 
    % Precompute inner temporal and spatial term & factors according to chosen model version
    ModelFuncs = set_modelFunctions( Inputs, ModelFuncs, sqrtParamETAS );
                                           
    % LL-function: 
    % First summand (Sum of logarithms of R0 at target events)
    idxFlag      = find(Catalog.flag > 0);
    [R0, dR0] = estimate_eventRates( Catalog, ...
                                        Inputs, ...
                                        sqrtParamETAS, ...
                                        ModelFuncs, ...
                                        idxFlag, ...
                                        -1, ...
                                        -1, -1, -1, -1, ...
                                        true, ...
                                        'time-space rates' ); 

    summand1_logL = sum( log(R0) .* (R0 > 1.0e-25) - 100 * (R0 <= 1.0e-25) );

    % Second summand (Integral of R0 over target space and time window)   
    [ LambdaIntegr, dLambdaIntegr, anisoPeaks ] = estimate_loglikIntegral( inclGrad, ...
                                                                           Catalog, ...
                                                                           Inputs, ...
                                                                           ModelFuncs, ...
                                                                           sqrtParamETAS, ...
                                                                           SpatData );

    summand2_logL           = sum( LambdaIntegr ) + sqrtParamETAS(1)^2 * aggregBackgrRate;

    % (-1) times sum of both 
    fval_logL               = - (summand1_logL - summand2_logL); 
    if ~isreal(fval_logL)
        test=-1;
    end

    % Calculate gradients
    if inclGrad

        % First term
        summand1_grad   = sum( dR0 ./ R0, 2 );

        % Second term
        summand2_grad   = sum( dLambdaIntegr, 2 ); 
        if ~Inputs.ModelSettings.fixedParamETAS(1)
            summand2_grad(1)= aggregBackgrRate * sqrtParamETAS(1) * 2;
        end

        % (-1) times sum of both 
        grad_logL       = - (summand1_grad - summand2_grad);

        % Display
        %mu_backgrPerDay = convert_mu2backgrPerDay( Catalog, TargetWindow, SpatData, sqrtParamETAS(1)^2 );
        mu_backgrPerDay = sqrtParamETAS(1)^2;

        % Display
        disp(['LL = ', num2str( round(-fval_logL, 3) )])
        disp(['LL_2 = ', num2str(summand2_logL), ' (', num2str(sum(Catalog.flag>0 & ~Catalog.isDupl)), ' target events)'])
%         round2dec = max(-floor(log10(abs(fval_logL))) + 2, 0);
%         fprintf(['Function Value = %0.', num2str(round2dec),'f\n'], fval_logL)

        print_iterationResults( sqrtParamETAS.^2, -grad_logL, Inputs.ModelSettings.fixedParamETAS )
        test=-1;

    else
        grad_logL = 'dummy';
        
    end
    
end
    