function stepLength = optimize_lineSearch( Catalog, ...
                                            Inputs, ...
                                            ModelFuncs, ...
                                            paramOld, ...
                                            searchDirection, ...
                                            loglik, ...
                                            aggregBackgrRate, ...
                                            stepLength, ...
                                            SpatData, ...
                                            TempData, ...
                                            writeLog )
    
    nFunctionCalls = 0;
    
    const2 = 1e-16;
    if stepLength <= 1e-30
        stepLength = 0.1;
    end
    
    searchDirecNorm = sqrt( sum(searchDirection.*searchDirection) );
    if searchDirecNorm > 1
        stepLength = stepLength/searchDirecNorm;
    end
    
    stepLength1 = 0;
    stepLength2 = stepLength;
    loglik1     = loglik;
    
    paramNew            = paramOld + stepLength2 * searchDirection;
    loglik2         = estimate_loglikSwitch(0, Catalog, Inputs, ModelFuncs, paramNew, aggregBackgrRate, SpatData, TempData, false);
    nFunctionCalls  = nFunctionCalls + 1;
    
    if abs(loglik2-loglik1) > 0.1*abs(loglik1)
        stepLength     = 0.05;
        if searchDirecNorm > 1
            stepLength = stepLength/searchDirecNorm;
        end
        stepLength2 = stepLength;
        paramNew        = paramOld + stepLength2 * searchDirection;
        loglik2     = estimate_loglikSwitch(0, Catalog, Inputs, ModelFuncs, paramNew, aggregBackgrRate, SpatData, TempData, false);
        nFunctionCalls = nFunctionCalls + 1;
        disp(['Reset stepLength2=0.05 since loglik2=', num2str(round(loglik2,3)), ' and loglik1=', num2str(round(loglik1,3))])
    end
    
    while loglik2 <= loglik1
%         disp(['loglik2=', num2str(loglik2), ' < loglik1=', num2str(loglik1), '; stepLengthda=', num2str(2*stepLength2)])
        stepLength3 = stepLength2 * 2;
        paramNew = paramOld + stepLength3 * searchDirection;
        loglik3 = estimate_loglikSwitch(0, Catalog, Inputs, ModelFuncs, paramNew, aggregBackgrRate, SpatData, TempData, false);
        nFunctionCalls = nFunctionCalls + 1;
%         disp(['Go further in direction. Log-likelihood improves by ', num2str(loglik3-loglik2)])
        if loglik3 > loglik2
            break;
        end
        stepLength1 = stepLength2;
        stepLength2 = stepLength3;
        loglik1 = loglik2;
        loglik2 = loglik3;
    end
        
    while loglik2>loglik1
        stepLength3 = stepLength2;
%         disp(['loglik2=', num2str(loglik2), ' > loglik1=', num2str(loglik1), '; stepLengthda=', num2str(0.1*stepLength3)])
        loglik3 = loglik2;
        stepLength2 = stepLength3 * 0.1;
        if stepLength2 * searchDirecNorm < const2
            stepLength = 0;
            return; % return what???
        end
        paramNew = paramOld + stepLength2 * searchDirection;
        loglik2 = estimate_loglikSwitch(0, Catalog, Inputs, ModelFuncs, paramNew, aggregBackgrRate, SpatData, TempData, false);
        nFunctionCalls = nFunctionCalls + 1;
%         disp(['Go shorter in direction. Log-likelihood improves by ', num2str(loglik2-loglik1)])
    end
    
    a1 = (stepLength3-stepLength2) * loglik1;
    a2 = (stepLength1-stepLength3) * loglik2;
    a3 = (stepLength2-stepLength1) * loglik3;
    b2 = (a1 + a2 + a3) * 2;
    b1 = a1 * (stepLength3+stepLength2) + a2 * (stepLength1+stepLength3) + a3 * (stepLength2+stepLength1);
    if b2==0
        stepLength = stepLength2;
        return;
    else
        stepLength = b1/b2;
        paramNew = paramOld + stepLength * searchDirection;
        loglik = estimate_loglikSwitch(0, Catalog, Inputs, ModelFuncs, paramNew, aggregBackgrRate, SpatData, TempData, false);
        nFunctionCalls = nFunctionCalls + 1;
        if stepLength > stepLength2
            if loglik <= loglik2
              stepLength1 = stepLength2;
              stepLength2 = stepLength;
              loglik1 = loglik2;
              loglik2 = loglik;
              % goto stat200130;
            else
              stepLength3 = stepLength;
              loglik3 = loglik;
              % goto stat200130;
            end
        else
            if (loglik >= loglik2)
              stepLength1 = stepLength;
              loglik1 = loglik;
              % goto stat200130;
            else
              stepLength3 = stepLength2;
              stepLength2 = stepLength;
              loglik3 = loglik2;
              loglik2 = loglik;
              % goto stat200130;
            end
        end
    end
    
    a1 = (stepLength3 - stepLength2)*loglik1;
    a2 = (stepLength1 - stepLength3)*loglik2;
    a3 = (stepLength2 - stepLength1)*loglik3;
    b2 = (a1 + a2 + a3)*2;
    b1 = a1 * (stepLength3 + stepLength2) + a2 * (stepLength1 + stepLength3) + a3 * (stepLength2 + stepLength1);
    if (b2 == 0)
        stepLength = stepLength2;
        return;
    else
        stepLength = b1 /b2;
        paramNew = paramOld + stepLength*searchDirection;
        loglik = estimate_loglikSwitch(0, Catalog, Inputs, ModelFuncs, paramNew, aggregBackgrRate, SpatData, TempData, false);
        nFunctionCalls = nFunctionCalls + 1;
        if (loglik2 < loglik)
          stepLength = stepLength2;
          return;
        end
    end
    
    if writeLog
        disp([num2str(nFunctionCalls), ' function calls.'])
        disp(['line search along the specified direction with step size = ', num2str(stepLength)])
    end
    
    if any(isnan(searchDirection)) || isnan(stepLength)
        error('NaN estimates in line search!');
    end

end

