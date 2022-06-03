function nAftershocks = sample_numberOfAftershocks( paramETAS, ...
                                                    SynthEventSet, ...
                                                    idxTriggerEvents, ...
                                                    SpatKernel, ...
                                                    f_inn, ...
                                                    g_inn, ...
                                                    M_c )

    % Extract parameters
    [~, A, alpha, c, p, ~, ~, q] = extract_paramETASvalues( paramETAS, 'default' );
    mag                          = SynthEventSet.mag(idxTriggerEvents);
    
    %% Compute temporal and spatial scaling
    % Temporal scaling: integral over Omori from 0 to T= min( tempExtent, nDays )
    tempExtent  = SynthEventSet.tempExtent(idxTriggerEvents);
    tempScaling   = 1/(1-p) * ( g_inn(tempExtent).^(1-p) - g_inn(0).^(1-p) );
    
    % Spatial scaling: integral over spatial kernel from 0 to spatExtent
    if SpatKernel.isNormalized
        spatScaling = 1; % ones(length(spatExtent),1);
    else
        spatExtent  = SynthEventSet.spatExtent(idxTriggerEvents);
        rupL        = SynthEventSet.rupL(idxTriggerEvents);
        wi          = SynthEventSet.wi(idxTriggerEvents);
        spatScaling = 1 - f_inn( spatExtent, spatExtent.^2, mag, rupL, wi ).^(1-q); 
    end
    
    %% Compute expected aftershock productivity
    expProductivity = A*exp(alpha*(mag-M_c)) .* tempScaling .* spatScaling .* SynthEventSet.evWeight(idxTriggerEvents);
    
    %% Sample number of aftershocks
    nAftershocks    = poissrnd( expProductivity );
    
end