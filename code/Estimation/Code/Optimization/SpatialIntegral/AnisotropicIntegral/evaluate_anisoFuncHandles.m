function [ iSample_integrAniso, ...
           iSample_integrIso, ...
           iSample_gradAniso_dD, ...
           iSample_gradIso_dD, ...
           iSample_gradAniso_dQ, ...
           iSample_gradIso_dQ ] = evaluate_anisoFuncHandles( distGrid_integr, ...
                                                            distGrid_grad, ...
                                                            F_integr, ...
                                                            F_integr_dD, ...
                                                            F_integr_dQ, ...
                                                            f_densAniso, ...
                                                            f_densAniso_dD, ...
                                                            f_densAniso_dQ, ...
                                                            f_densIso, ...
                                                            f_densIso_dD, ...
                                                            f_densIso_dQ, ...
                                                            isGradiants, ...
                                                            checkDetailed, ...
                                                            iMag, ...
                                                            iRupL, ...
                                                            iR )
                                                                
    %% Preallocations
    % Aniso integrals and gradiants
    aux_integrAniso     = zeros(length(distGrid_integr), 1);
    aux_gradAniso_dD    = zeros(length(distGrid_grad), 1);
    aux_gradAniso_dQ    = zeros(length(distGrid_grad), 1);
    
    % For checking purposes: iso integral and gradiants
    aux_integrIso   = zeros(length(distGrid_integr), 1);
    aux_gradIso_dD  = zeros(length(distGrid_grad), 1);    
    aux_gradIso_dQ  = zeros(length(distGrid_grad), 1);  
     
   %% Evaluate spatial integrals
    % Add copy of distance grid, with leading entry 0
    distGrid_integr0 = [0; distGrid_integr];
    % Loop over entries of iR
    for j = 2:length(distGrid_integr0)
        % Integrate anisotropic function between current distance grid points
        lowBound            = distGrid_integr0(j-1);
        uppBound            = distGrid_integr0(j);
%         r                   = distGrid_integr(j);            
        aux_integrAniso(j-1)  = integral(f_densAniso, lowBound, uppBound);
        
        % For checking purposes: Integrate isotropic function between current
        % distance grid points
        if checkDetailed
            aux_integrIso(j-1) = integral(f_densIso, lowBound, uppBound);
        end        
    end
    % Compute cumulative sum of integrals
    iSample_integrAniso = cumsum(aux_integrAniso);
    if checkDetailed
        checkSample_integrIso = cumsum(aux_integrIso);
    else
        checkSample_integrIso = -1;
    end
    
    % Compute isotropic integral as full analytical integral minus
    % anisotropic integral
    iSample_integrIso = F_integr(distGrid_integr) - iSample_integrAniso;

    %% Evaluate spatial gradiants
    if isGradiants
        
        distGrid_grad0 = [0; distGrid_grad];
        
        for j = 2:length(distGrid_grad0)
            % Integrate anisotropic function between current distance grid points
            lowBound            = distGrid_grad0(j-1);
            uppBound            = distGrid_grad0(j);
            aux_gradAniso_dD(j-1)  = integral(f_densAniso_dD, lowBound, uppBound);
            aux_gradAniso_dQ(j-1)  = integral(f_densAniso_dQ, lowBound, uppBound);
            
            % For checking purposes: Integrate isotropic gradiants between current distance grid
            % points
            if checkDetailed
                aux_gradIso_dD(j-1) = integral(f_densIso_dD, lowBound, uppBound);
                aux_gradIso_dQ(j-1) = integral(f_densIso_dQ, lowBound, uppBound);
            end
        end
        
        % Compute cumulative sum of integrals
        iSample_gradAniso_dD = cumsum(aux_gradAniso_dD);
        iSample_gradAniso_dQ = cumsum(aux_gradAniso_dQ);
        
        if checkDetailed
            checkSample_gradIso_dD = cumsum(aux_gradIso_dD);
            checkSample_gradIso_dQ = cumsum(aux_gradIso_dQ);
        else
            checkSample_gradIso_dD = -1;
            checkSample_gradIso_dQ = -1;
        end
        
        % Compute isotropic gradiants as full analytical gradiants minus
        % anisotropic gradiants
        iSample_gradIso_dD = F_integr_dD(distGrid_grad) - iSample_gradAniso_dD;
        iSample_gradIso_dQ = F_integr_dQ(distGrid_grad) - iSample_gradAniso_dQ;
        
    else
        
        iSample_gradIso_dD      = zeros(length(distGrid_grad), 1);
        iSample_gradIso_dQ      = zeros(length(distGrid_grad), 1);
        iSample_gradAniso_dD    = zeros(length(distGrid_grad), 1);
        iSample_gradAniso_dQ    = zeros(length(distGrid_grad), 1);
        checkSample_gradIso_dD  = -1;
        checkSample_gradIso_dQ  = -1;
        
    end

    %% Plausibility checks
    try
    check_anisoIntegrSamples( iSample_integrAniso, ...
                               iSample_integrIso, ...
                               iSample_gradAniso_dD, ...
                               iSample_gradIso_dD, ...
                               iSample_gradAniso_dQ, ...
                               iSample_gradIso_dQ, ...
                               checkSample_integrIso, ...
                               checkSample_gradIso_dD, ...
                               checkSample_gradIso_dQ, ...
                               iMag, ...
                               isGradiants, ...
                               checkDetailed )
    catch
        test=-1;
    end

    %% Checks I: Checking tools for more detailed investigation
    if checkDetailed && ( mod(iMag,1)==0 || iMag==8.7 )

        % Create plots
        plot_anisoIntegrEstim( iSample_integrAniso, ...
                               iSample_integrIso, ...
                               iSample_gradAniso_dD, ...
                               iSample_gradIso_dD, ...
                               iSample_gradAniso_dQ, ...
                               iSample_gradIso_dQ, ...
                               distGrid_integr, ...
                               distGrid_grad, ...
                               iMag, ...
                               iRupL, ...
                               iR )

    end
    
end