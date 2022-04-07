function print_iterationResults( paramETAS, gradients, fixedParamETAS )

    %% Parameter Conversions into Interpretable Units
    % Convert parameter D from degrees to km units
    paramETAS(6) = paramETAS(6) * 111^2;
    % Convert Tb from days to seconds unit
    paramETAS(9) = paramETAS(9) * 60*60*24;
    % Convert beta from exponential format to Gutenberg-Richter b
    paramETAS(10) = paramETAS(10) / log(10);
    
    %% Define format for display
    nLeadDigits_grad    = ceil(log10(max(abs(gradients))));
    strFormat_grad      = ['%', num2str( max(nLeadDigits_grad +5, 6) ), '.3f'];
    firstDigit_param    = ceil(log10(max(paramETAS)));   
    lastDigit_param     = floor( log10( min(paramETAS(paramETAS>0)) ) );
    strFormat_param     = ['%', num2str(firstDigit_param-lastDigit_param+3), '.', num2str(-lastDigit_param+2),'f'];
    
    %% Print iteration results
    % mu
    if ~fixedParamETAS(1)
        fprintf(['gradient[%d]  = ', strFormat_grad, '   mu    = ', strFormat_param, '\n'], ...
                1, gradients(1), paramETAS(1))
    end
    % A
    if ~fixedParamETAS(2)
        fprintf(['gradient[%d]  = ', strFormat_grad, '   A     = ', strFormat_param, '\n'], ...
                2, gradients(2), paramETAS(2))
    end
    % alpha
    if ~fixedParamETAS(3)
        fprintf(['gradient[%d]  = ', strFormat_grad, '   alpha = ', strFormat_param, ' (1/mag) \n'], ...
                3, gradients(3), paramETAS(3))
    end
    % c
    if ~fixedParamETAS(4)
        fprintf(['gradient[%d]  = ', strFormat_grad, '   c     = ', strFormat_param, ' (days) \n'], ...
                4, gradients(4), paramETAS(4))
    end
    % p
    if ~fixedParamETAS(5)
        fprintf(['gradient[%d]  = ', strFormat_grad, '   p     = ', strFormat_param, '\n'], ...
                5, gradients(5), paramETAS(5))
    end
    % D
    if ~fixedParamETAS(6)
        fprintf(['gradient[%d]  = ', strFormat_grad, '   D     = ', strFormat_param, ' (km^2) \n'], ...
                6, gradients(6), paramETAS(6))
    end
    % gamma
    if ~fixedParamETAS(7)
        fprintf(['gradient[%d]  = ', strFormat_grad, '   gamma = ', strFormat_param, ' (1/mag) \n'], ...
                7, gradients(7), paramETAS(7))
    end
    % q
    if ~fixedParamETAS(8)
        fprintf(['gradient[%d]  = ', strFormat_grad, '   q     = ', strFormat_param, '\n'], ...
                8, gradients(8), paramETAS(8))
    end
    % Tb
    if ~fixedParamETAS(9)
        fprintf(['gradient[%d]  = ', strFormat_grad, '   Tb    = ', strFormat_param, ' (sec) \n'], ...
                9, gradients(9), paramETAS(9))
    end
    % beta/ b
    if ~fixedParamETAS(10)
        fprintf(['gradient[%d] = ', strFormat_grad, '   b     = ', strFormat_param, '\n'], ...
                10, gradients(10), paramETAS(10))
    end

end

