function [ LambdaIntegral, ...
           dLambdaIntegral, ...
            anisoPeaks ] = estimate_loglikIntegral( inclGradients, ...
                                                    Catalog, ...
                                                    Inputs, ...
                                                    ModelFuncs, ...
                                                    sqrtParamETAS, ...
                                                    SpatData )

    %% Data preparation
    % Extract parameters (theta contains the sqrt of current parameter
    % values)
    [~,A,~,~,p,D,gamma,q] = extract_paramETASvalues( sqrtParamETAS, 'sqrt' );
    
    % Initialize vectors
    N               = length(Catalog.t);
    dLambdaIntegral = zeros(length(sqrtParamETAS),N);
    eventWeight     = Catalog.evWeight;
    
    %% Computation of integral and gradiants
    % Temporal integral and gradiants
    [ G, dG_c, dG_p ] = estimate_temporalIntegral( inclGradients, ...
                                                    Catalog, ...
                                                    p, ...
                                                    ModelFuncs.fTempK_inn, ...
                                                    1:N, ...
                                                    Inputs.TargetSettings.tWindow_days );
    
    % Spatial integral and gradiants
    [ F, dF_D, dF_gamma, dF_q, anisoPeaks ] = estimate_spatialIntegral( inclGradients, ...
                                                                          Catalog, ...
                                                                          Inputs, ...
                                                                          [D, gamma, q], ...
                                                                          ModelFuncs, ...
                                                                          SpatData );
    
    % Compute lambdaIntegral  
    Mc              = Inputs.TargetSettings.Mc;
    productivity    = ModelFuncs.fProductivity(Catalog.mag,Mc);
    prodXtime       = productivity .* G;
    LambdaIntegral  = prodXtime .* F;
    
    % Compute derivatives of lambdaIntegral
    if inclGradients

        prodXspace      = productivity .* F;
        fixedParamETAS  = Inputs.ModelSettings.fixedParamETAS;
        
        % Gradiants by productivity
        if ~fixedParamETAS(2)
            dLambdaIntegral(2,:)    = LambdaIntegral / A;
        end
        if ~fixedParamETAS(3)
            dLambdaIntegral(3,:)    = LambdaIntegral .* (Catalog.mag-Inputs.TargetSettings.Mc);
        end
        if ~fixedParamETAS(4)
            dLambdaIntegral(4,:)    = prodXspace .* dG_c;
        end
        if ~fixedParamETAS(5)
            dLambdaIntegral(5,:)    = prodXspace .* dG_p;
        end
        if ~fixedParamETAS(6)
            dLambdaIntegral(6,:)    = prodXtime .* dF_D;
        end
        if ~fixedParamETAS(7)
            dLambdaIntegral(7,:)    = prodXtime .* dF_gamma;
        end
        if ~fixedParamETAS(8)
            dLambdaIntegral(8,:)    = prodXtime .* dF_q;
        end
        
        % Account for inner derivative (since parameters are in squared form)
        dLambdaIntegral         = dLambdaIntegral * 2 .* sqrtParamETAS;
        
    end
    
    % Apply event weights
    LambdaIntegral  = LambdaIntegral .* eventWeight;
    dLambdaIntegral = dLambdaIntegral .* eventWeight';

end

