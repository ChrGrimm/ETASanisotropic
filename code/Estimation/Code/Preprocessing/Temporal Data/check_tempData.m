function check_tempData( Catalog, ...
                         Inputs, ...
                         ModelFuncs, ...
                         SpatData, ...
                         TempData, ...
                         sqrtParamETAS )

    %% Precomputations
    T2minusT1       = Inputs.TargetSettings.tWindow_days(2)-Inputs.TargetSettings.tWindow_days(1);
    
    % Precompute model functions according to chosen model version
    ModelFuncs = set_modelFunctions( Inputs, ModelFuncs, sqrtParamETAS );
                                       
    % Extract current parameter values
    [ ~, ~, ~, ~, ~, ...
      D, gamma, q, Tb ]      = extract_paramETASvalues( sqrtParamETAS, 'sqrt' );
  
    % Precompute spatial integral and gradiants
    computeGradients           	= ~all(Inputs.ModelSettings.fixedParamETAS(6:8));    
    [ F, dF_D, dF_gamma, dF_q ] = estimate_spatialIntegral( computeGradients, ...
                                                            Catalog, ...
                                                            Inputs, ...
                                                            [D, gamma, q], ...
                                                            ModelFuncs, ...
                                                            SpatData ); 
                                                                      
    %% Integral estimation
    nBackgrEventsPerDay = 0;
    tStart = tic;
    N0 = estimate_N0( Catalog, ...
                        Inputs, ...
                        TempData, ...
                        nBackgrEventsPerDay, ...
                        sqrtParamETAS, ...
                        ModelFuncs, ...
                        F, dF_D, dF_gamma, dF_q, ...
                        'with gradient' ); 
    duration = toc(tStart);
    disp(['Estimation of temporal integral takes roughly ', num2str(round(duration)), ' seconds.']) 
    
    % Compute second summand of LL-function
    expMinusN0  = exp(-N0);
    timeGrid    = TempData.timeGrid;
    LL_summand2 = T2minusT1/Tb - trapz( timeGrid, expMinusN0 ) / Tb;
    
%     dLL_summand2    = zeros(10,1);
%     dLL_summand2(1) = trapz( timeGrid, dN0(1,:) .* expMinusN0 ) / Tb;     % mu
%     dLL_summand2(2) = trapz( timeGrid, dN0(2,:) .* expMinusN0 ) / Tb;   % A
%     dLL_summand2(3) = trapz( timeGrid, dN0(3,:) .* expMinusN0 ) / Tb;   % alpha
%     dLL_summand2(4) = trapz( timeGrid, dN0(4,:) .* expMinusN0 ) / Tb;   % c
%     dLL_summand2(5) = trapz( timeGrid, dN0(5,:) .* expMinusN0 ) / Tb;   % p
%     dLL_summand2(6) = trapz( timeGrid, dN0(6,:) .* expMinusN0 ) / Tb;   % D
%     dLL_summand2(7) = trapz( timeGrid, dN0(7,:) .* expMinusN0 ) / Tb;   % gamma
%     dLL_summand2(8) = trapz( timeGrid, dN0(8,:) .* expMinusN0 ) / Tb;   % q
%     dLL_summand2(9) = ( -T2minusT1/Tb^2 + trapz( timeGrid, (Tb*dN0(9,:)+1) .* expMinusN0 ) / Tb^2 ) ...
%                         * 2 * sqrtParamETAS(9);   % Tb
%     dLL_summand2(10) = 0;                                               % beta
                        
    %% "Exact" integral (Matlab-internal integral estimation)
%     if exist([strCurrentDirectory, '\Exact_TempIntegral.mat'], 'file')
%         load('Exact_TempIntegral.mat')
%     else
    LL2_integExact = estimate_N0_integral( Catalog, ...
                                            Inputs, ...
                                            ModelFuncs, ...
                                            nBackgrEventsPerDay, ...
                                            sqrtParamETAS, ...
                                            F, dF_D, dF_gamma, dF_q, ...
                                            'with gradient' ); 
%         save Exact_TempIntegral.mat LL2_integExact dLL2_integExact LL2_exact
%     end
    
    disp('... Check temporal integral accuracy ...')
    deviationPercent = 100 * (round(LL_summand2/LL2_integExact, 3) - 1);
    disp(['Approximated integral varies from detailed estimation by ', num2str(deviationPercent), '%. Should be in the range of 0-2.5%.'])
    if deviationPercent > 2.5 && deviationPercent < 5
        warning('Approximation of temporal integral should deviate from detailed estimation by 0-2.5%. Consider modifying inputs in ''precompute_temporalData''.')
    elseif deviationPercent >= 5
        error('Approximation of temporal integral may not deviate from detailed estimation by more than 5%. Consider modifying inputs in ''precompute_temporalData''.')
    end
%     disp(['Ratio of Estimated Integral over "True" Integral = ', num2str(round(LL_summand2/LL2_integExact, 3))])                                                    

end