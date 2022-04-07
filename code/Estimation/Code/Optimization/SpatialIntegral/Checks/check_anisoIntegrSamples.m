function check_anisoIntegrSamples( iSample_integrAniso, ...
                                   iSample_integrIso, ...
                                   iSample_gradAniso_dD, ...
                                   iSample_gradIso_dD, ...
                                   iSample_gradAniso_dQ, ...
                                   iSample_gradIso_dQ, ...
                                   checkSample_integrIso, ...
                                   checkSample_gradIso_dD, ...
                                   checkSample_gradIso_dQ, ...
                                   iMag, ...
                                   isGradiant, ...
                                   checkDetailed )

    eps = 10^(-6);
    
    %% Check integral samples  
    % integral values must be positive
    if any( iSample_integrAniso < -eps ) || any( iSample_integrIso < -eps )
        error(['Integral sample for M=', num2str(iMag), ': iso and aniso interpolation values must be positive!'])
    % sum of iso and aniso integral must be smaller than 1
    elseif any( iSample_integrAniso + iSample_integrIso > 1 )
        error(['Integral sample for M=', num2str(iMag), ': Sum of iso and aniso interpolation values must be smaller than 1!'])
    % iso and aniso integral values must be strictly monotoneously increasing
    elseif any( diff(iSample_integrAniso) < -eps ) || any( diff(iSample_integrIso) < -eps )
        error(['Integral sample for M=', num2str(iMag), ': iso and aniso interpolation values must be strictly monotonically increasing!'])
    % iso and aniso integral value for r=10^(-30) must be below pre-defined eps
    elseif iSample_integrAniso(1) > eps || iSample_integrIso(1) > eps
        error(['Integral sample for M=', num2str(iMag), ': iso and aniso interpolation values for tiny distance must be smaller than eps!'])
%     % sum of iso and aniso integral for r=1000 must be larger than 1-eps
%     elseif (sample_F(end,1) + sample_F(end,2)) < 1-eps
%         error(['CDF sample_F ', num2str(iMag), ': iso and aniso interpolation values for huge distance must be larger than 1 - 10^(-4)!'])
    end
    
    if checkDetailed
        if any( abs(iSample_integrIso - checkSample_integrIso) > eps )
            error(['Integral sample for M=', num2str(iMag), ': iso sample is not consistent with check!'])
        end
    end
        
    %% Check gradiant samples
    if isGradiant
        
        %% Check gradiants with respect to D
        % aniso values must be smaller than 0
        if any( iSample_gradAniso_dD > 0 )
            error(['Gradiant sample by D for M=', num2str(iMag), ' by D: iso interpolation values must be negative!'])
        % Sum of iso and aniso values must be smaller than 0
        elseif any( iSample_gradAniso_dD + iSample_gradIso_dD > eps )
            error(['Gradiant sample by D for M=', num2str(iMag), ' by D: Sum of iso and aniso interpolation values must be negative!'])
        % absolute iso and aniso values for r=10^(-30) must be below pre-defined eps
        elseif abs(iSample_gradAniso_dD(1)) > eps || abs(iSample_gradIso_dD(1)) > eps
            error(['Gradiant sample by D for M=', num2str(iMag), ' by D: Absolute iso and aniso interpolation values for tiny distance must be smaller than eps=10^(-4)!'])
%         % sum of iso and aniso values for r=1000 must be larger than -eps
%         elseif sample_dF_D(end,1) + sample_dF_D(end,2) < -eps
%             error(['Gradiant sample by D for M=', num2str(iMag), ' by D: Sum of iso and aniso interpolation values for huge distance must be larger than -10^(-4)!'])
        end

        %% Check gradiants with respect to q        
        % aniso values must be larger than 0
        if any( iSample_gradAniso_dQ < 0 ) 
            error(['Gradiant sample by q for M=', num2str(iMag), ' by q: aniso interpolation values must be positive!'])
        % Sum of iso and aniso values must be larger than 0
        elseif any( iSample_gradAniso_dQ + iSample_gradIso_dQ < -eps )
            error(['Gradiant sample by q for M=', num2str(iMag), ' by q: Sum of iso and aniso interpolation values must be positive!'])
        % absolute iso and aniso values for r=10^(-30) must be below pre-defined eps
        elseif abs(iSample_gradAniso_dQ(1)) > eps || abs(iSample_gradIso_dQ(1)) > eps
            error(['Gradiant sample by q for M=', num2str(iMag), ' by q: Absolute iso and aniso interpolation values for tiny distance must be smaller than eps=10^(-4)!'])
%         % sum of iso and aniso values for r=1000 must be smaller than eps
%         elseif sample_dF_q(end,1) + sample_dF_q(end,2) > eps
%             error(['Gradiant sample by q for M=', num2str(iMag), ' by q: Sum of iso and aniso interpolation values for huge distance must be smaller than 10^(-4)!'])
        end
        
        if checkDetailed
            if any( abs(iSample_gradIso_dD - checkSample_gradIso_dD) > eps )
                error(['Gradiant sample by D for M=', num2str(iMag), ': iso sample is not consistent with check!'])
            end
            if any( abs(iSample_gradIso_dQ - checkSample_gradIso_dQ) > eps )
                error(['Gradiant sample by q for M=', num2str(iMag), ': iso sample is not consistent with check!'])
            end
        end
        
    end

end