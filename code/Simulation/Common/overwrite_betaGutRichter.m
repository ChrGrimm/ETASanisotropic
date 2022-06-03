function beta = overwrite_betaGutRichter( paramETAS, SimulSettings, ModelFuncs )

    %% If magnSampling is scalar, overwrite parameter beta by it
    if isscalar(SimulSettings.typeMagnitudeSampling)
        beta = log(10) * SimulSettings.typeMagnitudeSampling;
        disp(['Overwrite Gutenberg-Richter parameter by \beta = b*log(10) = ', num2str(beta)])
    else
        beta = paramETAS(10);
    end
    
    %% Compute and display branching ratio
    if isempty(SimulSettings.typeMagnitudeSampling) || isscalar(SimulSettings.typeMagnitudeSampling) 
        branchRatio = compute_branchingRatio( -1, -1, ...
                                              paramETAS(10), ...
                                              SimulSettings, ...
                                              ModelFuncs, ...
                                              'simulation' );
        disp(['Simulation is performed with branching ratio = ', num2str(branchRatio)])
    end

end

