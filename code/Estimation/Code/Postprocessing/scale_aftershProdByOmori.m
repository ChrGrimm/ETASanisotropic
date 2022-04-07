function [ fProductivity, ...
           fProductivity_scaled ] = scale_aftershProdByOmori( paramETAS, nDays )
    
    % Extract parameters
    A               = paramETAS(2);
    alpha           = paramETAS(3);
    c               = paramETAS(4);
    p               = paramETAS(5);
    
    % Compute Omori integral from 0 to nDays
    OmoriIntegral   = 1/(1-p) * ( (nDays+c).^(1-p) - c.^(1-p) );
    
    % Define function handle of scaled aftershock productivity 
    fProductivity           = @(mag,Mc) A * exp(alpha*(mag-Mc));
    fProductivity_scaled    = @(mag,Mc) OmoriIntegral * A * exp(alpha*(mag-Mc));                                          
                                               
end