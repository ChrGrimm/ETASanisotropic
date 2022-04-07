function [ fval_logL, ...
            grad_logL, ...
            LL_summand2, ...
            anisoPeaks] = estimate_loglikEtasIncomplete( inclGrad, ...
                                                         Catalog, ...
                                                         Inputs, ...
                                                         ModelFuncs, ...
                                                         sqrtParamETAS, ...
                                                         aggregBackgrRate, ...
                                                         SpatData, ...
                                                         TempData, ...
                                                         writeLog )
                                                     
% This function computes the LL-function (and optionally its derivatives) at the current iteration 
% for the ETAS-Incomplete model.    

    %% Extract data
    % Extract current parameter values
    [ ~, ~, ~, ~, ~, ...
      D, gamma, q, ...
      Tb, beta ]        = extract_paramETASvalues( sqrtParamETAS, 'sqrt' ); 
  
    % Extract column information from table Catalog
    idxFlag                 = find(Catalog.flag > 0);
    idxFlag_nonDuplicate    = idxFlag(~Catalog.isDupl(idxFlag));
    nTargetEvents           = length(idxFlag_nonDuplicate);
    mi                      = Catalog.mag(idxFlag_nonDuplicate);
    
    % Extract completeness magnitude
    Mc                      = Inputs.TargetSettings.Mc;
    
    % Extract length of target time interval
    T2minusT1               = Inputs.TargetSettings.tWindow_days(2)-Inputs.TargetSettings.tWindow_days(1);  
    
    %% Precomputations
    % Precompute model functions according to chosen model version
    ModelFuncs = set_modelFunctions( Inputs, ModelFuncs, sqrtParamETAS );
                                           
    % Precompute background event rate per day (in entire spatial target window)
    nBackgrEventsPerDay = aggregBackgrRate / T2minusT1;
    
    % Precompute spatial integral and gradiants
    computeGradients                        = inclGrad & ~all(Inputs.ModelSettings.fixedParamETAS(6:8));    
    [ F, dF_D, dF_gamma, dF_q, anisoPeaks ] = estimate_spatialIntegral( computeGradients, ...
                                                                        Catalog, ...
                                                                          Inputs, ...
                                                                          [D, gamma, q], ...
                                                                          ModelFuncs, ...
                                                                          SpatData );

    %% First Summand of Log-Likelihood Function
    % Compute "true" spatio-temporal event rate R0
    [R0, dR0] = estimate_eventRates( Catalog, ...
                                        Inputs, ...
                                        sqrtParamETAS, ...
                                        ModelFuncs, ...
                                        idxFlag, ...
                                        -1, ...
                                        -1, -1, -1, -1, ...
                                        true, ...
                                        'time-space rates' ); 
                                
    % Compute logarithm of "true" event rate R0                            
    logR0       = log(R0) .* (R0 > 1.0e-25) - 100 * (R0 <= 1.0e-25);
                                
    % Compute "true" temporal event rate N0 (over entire spatial target window)
    [N0, dN0]   = estimate_eventRates( Catalog, ...
                                        Inputs, ...
                                        sqrtParamETAS, ...
                                        ModelFuncs, ...
                                        idxFlag, ...
                                        nBackgrEventsPerDay, ...
                                        F, dF_D, dF_gamma, dF_q, ...
                                        true, ...
                                        'time rates' );
    
    % Compute first summand of LL-function
    magnTerm                = beta*(mi-Mc)';
    LL_summand1             = nTargetEvents * log(beta) + sum( logR0 - magnTerm - N0 .* exp(-magnTerm) );
    
    % Compute gradients of first summand of LL-function
    if inclGrad
        dLL_summand1            = zeros(10,1);
        dLL_summand1(1)         = sum( dR0(1,:)./R0 - dN0(1,:) .* exp(-magnTerm) );  % mu
        dLL_summand1(2)         = sum( dR0(2,:)./R0 - dN0(2,:) .* exp(-magnTerm) );  % A
        dLL_summand1(3)         = sum( dR0(3,:)./R0 - dN0(3,:) .* exp(-magnTerm) );  % alpha
        dLL_summand1(4)         = sum( dR0(4,:)./R0 - dN0(4,:) .* exp(-magnTerm) );  % c
        dLL_summand1(5)         = sum( dR0(5,:)./R0 - dN0(5,:) .* exp(-magnTerm) );  % p
        dLL_summand1(6)         = sum( dR0(6,:)./R0 - dN0(6,:) .* exp(-magnTerm) );  % D
        dLL_summand1(7)         = sum( dR0(7,:)./R0 - dN0(7,:) .* exp(-magnTerm) );  % gamma
        dLL_summand1(8)         = sum( dR0(8,:)./R0 - dN0(8,:) .* exp(-magnTerm) );  % q
        dLL_summand1(9)         = sum( - dN0(9,:) .* exp(-magnTerm) );               % Tb
        if ~Inputs.ModelSettings.fixedParamETAS(10)
            dLL_summand1(10)        = nTargetEvents/beta + sum( (mi-Mc)' .* ( N0.*exp(-magnTerm) - 1) );   % beta
        end
    end
    
    %% Second Summand of Log-Likelihood Function
    % Estimate "true" temporal event rates N0 at time grid (over entire spatial target window)
    [N0,dN0] = estimate_N0( Catalog, ...
                            Inputs, ...
                            TempData, ...
                            nBackgrEventsPerDay, ...
                            sqrtParamETAS, ...
                            ModelFuncs, ...
                            F, dF_D, dF_gamma, dF_q, ...
                            'with gradient' ); 

    % Compute second summand of LL-function
    expMinusN0  = exp(-N0);
    timeGrid    = TempData.timeGrid;
    LL_summand2 = T2minusT1/Tb - trapz( timeGrid, expMinusN0 ) / Tb; 
    
    if inclGrad
        dLL_summand2    = zeros(10,1);
        dLL_summand2(1) = trapz( timeGrid, dN0(1,:) .* expMinusN0 ) / Tb;     % mu
        dLL_summand2(2) = trapz( timeGrid, dN0(2,:) .* expMinusN0 ) / Tb;   % A
        dLL_summand2(3) = trapz( timeGrid, dN0(3,:) .* expMinusN0 ) / Tb;   % alpha
        dLL_summand2(4) = trapz( timeGrid, dN0(4,:) .* expMinusN0 ) / Tb;   % c
        dLL_summand2(5) = trapz( timeGrid, dN0(5,:) .* expMinusN0 ) / Tb;   % p
        dLL_summand2(6) = trapz( timeGrid, dN0(6,:) .* expMinusN0 ) / Tb;   % D
        dLL_summand2(7) = trapz( timeGrid, dN0(7,:) .* expMinusN0 ) / Tb;   % gamma
        dLL_summand2(8) = trapz( timeGrid, dN0(8,:) .* expMinusN0 ) / Tb;   % q
        if ~Inputs.ModelSettings.fixedParamETAS(9)
            dLL_summand2(9) = ( -T2minusT1/Tb^2 + trapz( timeGrid, (Tb*dN0(9,:)+1) .* expMinusN0 ) / Tb^2 ) ...
                                * 2 * sqrtParamETAS(9);   % Tb
        end
        dLL_summand2(10) = 0;                                               % beta
    end       
    
    %% Put Together Summands of Log-Likelihood Function
    fval_logL               = - (LL_summand1 - LL_summand2); 
              
    if inclGrad
        grad_logL           = - (dLL_summand1 - dLL_summand2);
        
        if writeLog
            % Print current parameter estimates and gradients
    %         round2dec       = 3; max(-floor(log10(abs(fval_logL))) + 2, 0);
            disp(['LL = ', num2str( round(-fval_logL, 3) )])
            disp(['LL_2 = ', num2str(LL_summand2), ' (', num2str(sum(Catalog.flag>0 & ~Catalog.isDupl)), ' target events)'])
        
            print_iterationResults( sqrtParamETAS.^2, -grad_logL, Inputs.ModelSettings.fixedParamETAS )
        end
        
    else
        % dummy to assign output value
        grad_logL = 'dummy';
        
    end
        
end