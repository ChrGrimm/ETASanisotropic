function IterationLog = update_iterationLog( IterationLog, ...
                                                loglik, ...
                                                loglik_gradient, ...
                                                integratedRate, ...
                                                paramETAS, ...
                                                iIter )

    
    %% Update diagnostic data
    IterationLog.loglik_fval(iIter,1)   = -loglik;
    IterationLog.loglik_grad(iIter,:)   = loglik_gradient;
    nFreeParameters                     = sum( ~(paramETAS==0) );
    IterationLog.aic(iIter,1)           = 2 * (loglik + nFreeParameters);
    IterationLog.paramETAS(iIter,:)     = paramETAS';
    IterationLog.integratedRate(iIter,1)= integratedRate;

    disp(['loglikelihood = ', num2str(-loglik), ' ;    AIC = ', num2str(2*(loglik+8))]);

end

%% Old code pieces
%     IterationLog.avcov                = 1/4 * diag(1./sqrtParamETAS) * hessianMatrix * diag(1./sqrtParamETAS);
%     IterationLog.asd(iIter,:)         = sqrt(diag(IterationLog.avcov))';
%         round2dec_fv_old    = max(-floor(log10(abs(IterationLog.loglik_fval(iIter-1)))) + 2, 0);
%         round2dec_fv_new    = max(-floor(log10(abs(loglik))) + 2, 0);
%         round2dec_p_old     = max(-floor(log10(IterationLog.paramETAS(iIter-1,:))) + 2, 0);
%         round2dec_p_new     = max(-floor(log10(sqrtParamETAS.^2)) + 2, 0);
%         
%         round( IterationLog.loglik_fval(iIter-1), round2dec_fv_old ) == round( -loglik, round2dec_fv_new ) ...
%                 && all( round2dec_p_old == round2dec_p_new' ) ...
%                 && all( round( IterationLog.paramETAS(iIter-1,:) .* 10.^round2dec_p_old ) ...
%                         == round( (sqrtParamETAS.^2 .* 10.^round2dec_p_new)' ) )
