function ModelFuncs = set_modelFunctions( Inputs, ...
                                          ModelFuncs, ...
                                          sqrtParamETAS )
%
% This function prepares function handles for the inner parts of the
% temporal and spatial triggering distributions as well as the scaling
% from both temporal and spatial distribution and productivity term of the 
% ETAS model for various model designs
%
% Call: [ g_inn, f_inn, factor ] = set_modelFunctions( modelDesign, sqrtParamETAS, magn )
%
% INPUTS:
%
% sqrtParamETAS:  numeric vector (8x1)    square roots of estimated ETAS parameters
% modelDesign:  string variable         controlling the chosen model design
% Mc:           scalar                  (lower) magnitude threshold

    %% Extract parameter values
    [~, A, alpha, c, p, D, gamma, q] = extract_paramETASvalues( sqrtParamETAS, 'sqrt' );
    Mc = Inputs.TargetSettings.Mc;
    
    %% Spatial kernel
    % Define spatial functions
    if strcmp(Inputs.ModelSettings.spaceModel, 'none')
        ModelFuncs.fSpatK_inn       = @(r, r2, mag, rupExtent) 1;
        ModelFuncs.fSpatK_factor    = @(R, mag, rupExtent) 1;
        ModelFuncs.fSpatK_width     = @(mag) 1;
        
    else
        ModelFuncs.fSpatK_inn       = @(r, r2, mag, rupExtent) 1 + (2*rupExtent.*r + pi*r2) ./ (D*exp(gamma*(mag-Mc)));
        ModelFuncs.fSpatK_factor    = @(R, mag, rupExtent) (q-1) ./ ( D*exp(gamma*(mag-Mc)) .* (1 - ModelFuncs.fSpatK_inn(R,R.^2,mag,rupExtent).^(1-q)) );
        ModelFuncs.fSpatK_width     = @(mag) D*exp(gamma*(mag-Mc));

    end
    
    %% Temporal design
    % Define temporal functions
    ModelFuncs.fTempK_inn       = @(t)  t + c;
    ModelFuncs.fTempK_factor    = @(T)  1;

    %% Productivity and factor
    ModelFuncs.fProductivity        = @(mag, Mc) A * exp(alpha*(mag-Mc));
    nDays                           = Inputs.TimeSettings.restr_days;
    OmoriIntegral                   = 1/(1-p) * ( (nDays+c).^(1-p) - c.^(1-p) );
    ModelFuncs.fProductivity_scaled = @(mag, Mc) OmoriIntegral * A * exp(alpha*(mag-Mc));
            
end

%% Old code pieces
%     if strcmp( Inputs.ModelSettings.timeModel, 'Omori' )
        
%     else
%         ModelFuncs.fTempK_inn       = 1 + t./c;
%         ModelFuncs.fTempK_factor    = (p-1)./c ./ (1 - g_inn(T).^(1-p));
%     end
        
%     if ismember(SpatKernel.type, {'default', 'dip'})
%         f_inn           = @(r,r2,mi,rupLi)    1 + (2*rupLi.*r + pi*r2) ./ (D*exp(gamma*(mi-Mc)));
%         sigma           = @(mag)                 D*exp(gamma*(mag-Mc));
%         % spatial factor
%         spatialFactor   = @(R,mi,rupLi)   (q-1) ./ ( D*exp(gamma*(mi-Mc)) .* (1 - f_inn(R,R.^2,mi,rupLi).^(1-q)) );
%         
%     elseif strcmp(SpatKernel.type, 'none')
%         f_inn           = @(r,r2,mi,rupLi)   1;
%         sigma           = @(mag)                1;
%         spatialFactor   = @(R,mi,rupLi)      1;
%         
%     else
%         error('Unknown version of spatial kernel')