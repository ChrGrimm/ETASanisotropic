function [LL2_integ, ...
            dLL2_integ, ...
            LL2_values] = estimate_N0_integral( Catalog, ...
                                                Inputs, ...
                                                ModelFuncs, ...
                                                nBackgrPerDay, ...
                                                sqrtParamETAS, ...
                                                F, dF_D, dF_gamma, dF_q, ...
                                                strMode )

    computeGradients    = strcmp(strMode, 'with gradient');
    
    % Extract ETAS parameters (they are in square root format  for numerical reasons)
    [mu, A, ~, c, p, ~, ~, ~, Tb] = extract_paramETASvalues( sqrtParamETAS, 'sqrt' );
    
    % Extract event-specific information in columns of table Catalog
    t           = Catalog.t;
    mag         = Catalog.mag;
    evWeight    = Catalog.evWeight;
    tempRestr   = Catalog.tempRestr;
    
    Mc          = Inputs.TargetSettings.Mc;
    
    factor          = ModelFuncs.fProductivity(mag, Mc) .* evWeight .* F;
%     factor_dD       = productivity .* evWeight .* dF_D;
%     factor_dGamma   = productivity .* evWeight .* dF_gamma;
%     factor_dQ       = productivity .* evWeight .* dF_q;
    
    idxTargetTime   = find(Catalog.flag>=0);
    timeGrid        = [ Inputs.TargetSettings.tWindow_days(1,1); ...
                        t(idxTargetTime); ...
                        Inputs.TargetSettings.tWindow_days(end,2) ];
    LL2_integ       = zeros(length(idxTargetTime)-1, 1);
    dLL2_integ      = zeros(length(idxTargetTime)-1, 10);
    LL2_values      = zeros(length(idxTargetTime)-1, 1);
    
    tic
    for jCounter = 1:length(idxTargetTime)-1
        
        tStart  = timeGrid(jCounter);
        tEnd    = timeGrid(jCounter+1)-10^-10;
        
        isContribEv = t <= tStart & t >= tStart-tempRestr; % & t >= tStart-0.1;
        tContr      = t(isContribEv);
        factorContr = factor(isContribEv);
        if ~isempty(tContr)
            if length(tContr)==1
                tContr      = [tContr; tContr];
                factorContr = [factorContr; 0];
            end
            funcN0_trig = @(tt) factorContr .* (tt-tContr+c).^(-p);
            funcN0      = @(tt) Tb * ( mu * nBackgrPerDay + sum( funcN0_trig(tt) ) );
            funcLL2      = @(tt) (1 - exp(- funcN0(tt) )) / Tb;
%             funcLL2_mu   = @(tt) exp(- funcN0(tt) ) * nBackgrPerDay * 2 * sqrtParamETAS(1);
%             funcLL2_A    = @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt)/A ) * 2 * sqrtParamETAS(2);
%             funcLL2_alpha= @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt) .* (mag(isContribEv)-Mc) ) * 2 * sqrtParamETAS(3);
%             funcLL2_c    = @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt) .* (-p./(tt-tContr+c)) ) * 2 * sqrtParamETAS(4);
%             funcLL2_p    = @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt) .* (-log(tt-tContr+c)) ) * 2 * sqrtParamETAS(5);
%             funcLL2_Tb   = @(tt) ((-1/Tb^2) + exp(- funcN0(tt) ) / Tb^2 .* ( Tb * (mu * nBackgrPerDay + sum( funcN0_trig(tt) )) + 1)) * 2 * sqrtParamETAS(9);

            LL2_integ(jCounter)     = integral(funcLL2, tStart, tEnd, 'RelTol', 1e-03);
%             dLL2_integ(jCounter,1)  = integral(funcLL2_mu, tStart, tEnd, 'RelTol', 1e-03);
%             dLL2_integ(jCounter,2)  = integral(funcLL2_A, tStart, tEnd, 'RelTol', 1e-03);
%             dLL2_integ(jCounter,3)  = integral(funcLL2_alpha, tStart, tEnd, 'RelTol', 1e-03);
%             dLL2_integ(jCounter,4)  = integral(funcLL2_c, tStart, tEnd, 'RelTol', 1e-03);
%             dLL2_integ(jCounter,5)  = integral(funcLL2_p, tStart, tEnd, 'RelTol', 1e-03);
%             dLL2_integ(jCounter,9)  = integral(funcLL2_Tb, tStart, tEnd, 'RelTol', 1e-03);

        end
    end
    toc
    
% %     tic
%     for jCounter = 1:length(timeGrid)-1
%         
%         tStart  = timeGrid(jCounter);
%         tEnd    = timeGrid(jCounter+1)-10^-10;
%         
%         tStep       = 100; %nthroot(160./(p*(p+1)*(tEnd-tStart)), -p-2) -c; % max(tEnd-tStart, 0.1);
%         isContribEv = t <= tStart & (t >= tStart-tStep | mag > 5);
%         isBefore    = t < tStart-tStep & mag <= 5;
%         tContr      = t(isContribEv);
%         tBefore     = t(isBefore);
% 
%         factorContr = factor(isContribEv);
%         factorBefore= factor(isBefore);
%         if isempty(tContr)
%             N0_2(jCounter) = 0;
%         else
%             if length(tContr)==1
%                 tContr      = [tContr; tContr];
%                 factorContr = [factorContr; 0];
%             end
%             N0_before_start = sum(factorBefore .* (tStart-tBefore+c).^(-p));
%             N0_before_end   = sum(factorBefore .* (tEnd-tBefore+c).^(-p));
%             m_N0_before     = (N0_before_end - N0_before_start) / (tEnd-tStart);
%             
%             funcN0_trigContr    = @(tt) factorContr .* (tt-tContr+c).^(-p);
%             funcN0_trigBefore   = @(tt) N0_before_start + m_N0_before * (tt-tStart);
%             funcN0_trig         = @(tt) sum(funcN0_trigContr(tt)) + funcN0_trigBefore(tt);
%             funcN0              = @(tt) Tb * ( mu * nBackgrPerDay + funcN0_trig(tt) );
%             funcLL2             = @(tt) ( 1 - exp(-funcN0(tt)) ) / Tb;
% %             funcLL2_mu          = @(tt) exp(- funcN0(tt) ) * nBackgrPerDay * 2 * sqrtParamETAS(1);
% %             funcLL2_A           = @(tt) exp(- funcN0(tt) ) .* funcN0_trig(tt) / A * 2 * sqrtParamETAS(2);
% %             funcLL2_alpha       = @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt) .* (mag(isContribEv)-Mc) ) * 2 * sqrtParamETAS(3);
% %             funcLL2_c           = @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt) .* (-p./(tt-tContr+c)) ) * 2 * sqrtParamETAS(4);
% %             funcLL2_p           = @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt) .* (-log(tt-tContr+c)) ) * 2 * sqrtParamETAS(5);
% %             funcLL2_Tb          = @(tt) ((-1/Tb^2) + exp(- funcN0(tt) ) / Tb^2 .* ( Tb * (mu * nBackgrPerDay + sum( funcN0_trig(tt) )) + 1)) * 2 * sqrtParamETAS(9);
% 
%             LL2_integ(jCounter) = integral(funcLL2, tStart, tEnd, 'AbsTol', 1e-03);
%             LL2_values(jCounter)= funcLL2( tStart );
%         end
%     end
%     toc
    
    LL2_integ = sum(LL2_integ);
    dLL2_integ = sum(dLL2_integ);
    
%     idxTimePoint    = TempData.idxTimePoint;
%     idxContrEvents  = TempData.idxContrEvents;
%     
    
end
