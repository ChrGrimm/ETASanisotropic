function ModelSettings = preprocess_modelSettings( ModelSettings )

    %% Extract and convert initial parameter guesses
    ModelSettings.iniParamETAS = [ ModelSettings.mu_ini, ModelSettings.A_ini, ...
                                   ModelSettings.alpha_ini, ModelSettings.c_ini, ...
                                   ModelSettings.p_ini, ModelSettings.D_ini, ...
                                   ModelSettings.gamma_ini, ModelSettings.q_ini, ...
                                   ModelSettings.Tb_ini, ModelSettings.b_ini ]';
    % Convert D to unit degrees
    ModelSettings.iniParamETAS(6) = ModelSettings.iniParamETAS(6)/111^2;
    % Convert Tb to unit days
    ModelSettings.iniParamETAS(9) = ModelSettings.iniParamETAS(9)/(60*60*24);
    % Convert b-value to exponential basis
    ModelSettings.iniParamETAS(10) = ModelSettings.iniParamETAS(10) * log(10);
    
    if ~isfield(ModelSettings, 'iniBackgrProb')
        ModelSettings.iniBackgrProb = 1;
    end
    
    %% Modify model settings
    ModelSettings.fixedParamETAS = false(1,10);
    % Set unused parameters to -1
    if strcmp( ModelSettings.spaceModel, 'none' )
%         ModelSettings.iniParamETAS(6:8)  = NaN;
        ModelSettings.fixedParamETAS(6:8)= true;
    end
    if ~strcmp( ModelSettings.modelType, 'ETAS-Incomplete' )
%         ModelSettings.iniParamETAS(9:10) = NaN;
        ModelSettings.fixedParamETAS(9:10)= true;
    end
    % Set p=1.000001 if p=1 for numerical reasons
    if ModelSettings.iniParamETAS(5)==1
        ModelSettings.iniParamETAS(5) = ModelSettings.iniParamETAS(5) + 10^-6;
    end
    % isParameterFixed must be row vector
    if size(ModelSettings.fixedParamETAS,1) > 1
        ModelSettings.fixedParamETAS = ModelSettings.fixedParamETAS';
    end
    
    %% Check plausibility of inputs
    % Fulfills constraints
    if any( ~(length(ModelSettings.iniParamETAS)==10) | ~isnumeric(ModelSettings.iniParamETAS) )
        error('iniParamETAS must be a numeric vector of length 10 with numeric entries.')
    end
    
end