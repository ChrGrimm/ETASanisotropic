function [ R0, ...
            dR0, ...
            R0_matrix, ...
            N0, ...
            dN0 ] = initialize_eventRateOutputs( nFlagEvents, ...
                                                 nParamETAS, ...
                                                 nEventsTotal, ...
                                                 computeGradients, ...
                                                 computeTriggerMatrix )
                                                         
    % Event rate vector is needed in all modes
    R0 = zeros( 1, nFlagEvents);
    N0 = zeros( 1, nFlagEvents);             
    % Gradient matrix 
    if computeGradients
        dR0 = zeros( nParamETAS, nFlagEvents);
        dN0 = zeros( nParamETAS, nFlagEvents);
    else
        dR0 = -1; % dummy
        dN0 = -1;
    end
    if computeTriggerMatrix
        R0_matrix  = zeros( nEventsTotal, nFlagEvents);
    else
        R0_matrix  = -1; % dummy
    end                                                     
                                                         
end