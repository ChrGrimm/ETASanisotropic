function [ F_integral, ...
           F_integral_dD, ...
           F_integral_dQ, ...
           f_densityAniso, ...
           f_densityAniso_dD, ...
           f_densityAniso_dQ, ...
           f_densityIso, ...
           f_densityIso_dD, ...
           f_densityIso_dQ ] = define_anisoFuncHandles( fSpatK_inn, ...
                                                        fSpatK_factor, ...
                                                        spatParamETAS, ...
                                                        R, ...
                                                        mag, ...
                                                        rupExtent, ...
                                                        isGradiants, ...
                                                        checkDetailed )
       
    D       = spatParamETAS(1);
    q       = spatParamETAS(3);
    
    %% Preallocate dummies for undefined functions
    if ~isGradiants 
        F_integral_dD       = -1;
        F_integral_dQ       = -1;
        f_densityAniso_dD   = -1;
        f_densityAniso_dQ   = -1;
        f_densityIso_dD     = -1;
        f_densityIso_dQ     = -1;
    end
    if ~checkDetailed
        f_densityIso        = -1;
        f_densityIso_dD     = -1;
        f_densityIso_dQ     = -1;
    end
    
    %% Precompute normalization factor and its contribution to gradiants
    f_inn_R             = fSpatK_inn(R, R.^2, mag, rupExtent);
    f_inn_R_to1minusQ   = f_inn_R.^(1-q);
    normalizFactor      = 1 - f_inn_R_to1minusQ;
%         gradIntegral       = [ 
    normaliz_integr_dD  = f_inn_R.^(-q) .* (f_inn_R - 1);
    normaliz_integr_dQ  = log(f_inn_R) .* f_inn_R_to1minusQ;
    normaliz_density_dD = (1-q) * normaliz_integr_dD / normalizFactor;
    normaliz_density_dQ = normaliz_integr_dQ / normalizFactor;
    
    %% Analytical integral and density function
    F_integral      = @(r) (1 - fSpatK_inn(r, r.^2, mag, rupExtent).^(1-q)) / normalizFactor;
    f_densityAniso  = @(r) fSpatK_factor(R, mag, rupExtent) * (2*rupExtent) .* fSpatK_inn(r, r.^2, mag, rupExtent).^(-q);
    if checkDetailed
       f_densityIso = @(r) fSpatK_factor(R, mag, rupExtent) * (2*pi*r) .* fSpatK_inn(r, r.^2, mag, rupExtent).^(-q);
    end

    %% Spatial gradiant functions
    if isGradiants
        % Normalization factor squared
        normalizSqu = normalizFactor^2; 
        
        % Derivative by D
        % analytical integral
        F_integral_dD       = @(r) (1-q)/(D*normalizSqu) ...
                                .* ( fSpatK_inn(r, r.^2, mag, rupExtent).^(1-q) .* (1-fSpatK_inn(r, r.^2, mag, rupExtent).^(-1)) * normalizFactor ...
                                     - (1-fSpatK_inn(r, r.^2, mag, rupExtent).^(1-q)) * normaliz_integr_dD );
        % anisotropic density
        f_densityAniso_dD   = @(r) f_densityAniso(r) * 1/D .* ...
                                    (q * (1-fSpatK_inn(r, r.^2, mag, rupExtent).^(-1)) - 1 - normaliz_density_dD);
        % isotropic density
        if checkDetailed
            f_densityIso_dD = @(r) f_densityIso(r) * 1/D .* ...
                                    (q * (1-fSpatK_inn(r, r.^2, mag, rupExtent).^(-1)) - 1 - normaliz_density_dD); 
        end    

        % Derivative by q
        % analytical integral
        F_integral_dQ       = @(r) 1/normalizSqu ...
                                .* ( log(fSpatK_inn(r, r.^2, mag, rupExtent)) .* fSpatK_inn(r, r.^2, mag, rupExtent).^(1-q) * normalizFactor ...
                                     - (1-fSpatK_inn(r, r.^2, mag, rupExtent).^(1-q)) * normaliz_integr_dQ );
        % anisotropic density
        f_densityAniso_dQ   = @(r) f_densityAniso(r) .* ...
                                    (1/(q-1) - log(fSpatK_inn(r, r.^2, mag, rupExtent)) - normaliz_density_dQ);
        % isotropic density
        if checkDetailed
            f_densityIso_dQ = @(r) f_densityIso(r) .* ...
                                    (1/(q-1) - log(fSpatK_inn(r, r.^2, mag, rupExtent)) - normaliz_density_dQ);
        end
            
    end

end

%% Old code pieces
% else
%     if isWiderThanMinWidth
%         %% Derivative by D
%         % analytical integral
%         F_integral_dD       = @(r) (1-q)/D ...
%                                 * fSpatK_inn(r, r.^2, mag, rupExtent).^(1-q) .* (1-fSpatK_inn(r, r.^2, mag, rupExtent).^(-1));
%         % anisotropic density
%         f_densityAniso_dD   = @(r) f_densityAniso(r) * 1/D .* ...
%                                     (q * (1-fSpatK_inn(r, r.^2, mag, rupExtent).^(-1)) - 1);
%         % isotropic density
%         if checkDetailed
%             f_densityIso_dD = @(r) f_densityIso(r) * 1/D .* ...
%                                     (q * (1-fSpatK_inn(r, r.^2, mag, rupExtent).^(-1)) - 1);  
%         end
%     end
% 
%     %% Derivative by q
%     % analytical integral
%     F_integral_dQ       = @(r) log(fSpatK_inn(r, r.^2, mag, rupExtent)) .* fSpatK_inn(r, r.^2, mag, rupExtent).^(1-q);
%     % anisotropic density
%     f_densityAniso_dQ   = @(r) f_densityAniso(r) .* ...
%                                 (1/(q-1) - log(fSpatK_inn(r, r.^2, mag, rupExtent)));
%     % isotropic density
%     if checkDetailed
%         f_densityIso_dQ = @(r) f_densityIso(r) .* ...
%                                 (1/(q-1) - log(fSpatK_inn(r, r.^2, mag, rupExtent)));
%     end
% 
% end


