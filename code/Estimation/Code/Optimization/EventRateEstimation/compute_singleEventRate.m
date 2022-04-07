function [ eventRate, ...
           eventRate_grad, ...
           triggRateMatrix ] = compute_singleEventRate( iTargetEvent, ...
                                                         t, mag, tempRestr, evWeight, isInExtent, dist, dist2, rupExtent, backgrRate, ...
                                                         fSpatK_inn, fTempK_inn, ...
                                                         mu, A, p, D, q, Tb, ...
                                                         Mc, nBackgrPerDay, spaceModel, ...
                                                         norm_D, norm_q, ...
                                                         factor, factor_dD, factor_dGamma, factor_dQ, ...
                                                         computeGradients, fixedParamETAS, ...
                                                         strModus )
              
    eventRate_grad = zeros(10,1);
    
    %% Compute triggered rate
    % Precomputations for time-space version
    if ismember(strModus, {'time-space rates', 'decluster'})
        % Extract contributing events
        isContributing  = isInExtent{iTargetEvent};
        % Evaluate spatial kernel
        fSpatK_inn_ij   = fSpatK_inn( dist{iTargetEvent}, dist2{iTargetEvent}, ...
                                      mag(isContributing), rupExtent(isContributing) );
        

    % Precomputations for time version
    elseif strcmp( strModus, 'time rates')
        % Determine contributing events
        isContributing  = t < t(iTargetEvent) & ...
                          t >= t(iTargetEvent) - tempRestr; % & t(iTargetEvent)-t <= 365;
        fSpatK_inn_ij   = 1;

    end
    
    % Compute weights
    weights         = evWeight(isContributing) * evWeight(iTargetEvent);
    
    % Evaluate time kernel
    fTempK_inn_ij           = fTempK_inn(t(iTargetEvent)-t(isContributing));
    fTempK_inn_ij_toMinusP  = fTempK_inn_ij.^(-p);

    % Compute and sum up event rates
    triggRates_ij  = weights .* factor(isContributing) .* fTempK_inn_ij_toMinusP .* fSpatK_inn_ij.^(-q);
    triggRate      = sum( triggRates_ij );

    % Fill trigger rate matrix
    triggRateMatrix                 = zeros(size(t));
    if ~isempty(triggRates_ij) 
        triggRateMatrix(isContributing) = triggRates_ij;
    end
    
    
    %% Compute event rate
    if ismember(strModus, {'time-space rates', 'decluster'})
        eventRate = triggRate + mu * backgrRate(iTargetEvent);
        
    elseif strcmp( strModus, 'time rates')
        eventRate = Tb * (mu * nBackgrPerDay + triggRate);
    end

    
    %% Gradients
    if computeGradients
        % by mu
        if ismember(strModus, {'time-space rates', 'decluster'})
            eventRate_grad(1)   = backgrRate(iTargetEvent);
            Tb                  = 1;
        elseif strcmp(strModus, 'time rates')
            if strcmp( spaceModel, 'none' )
                eventRate_grad(1) = Tb;
            else
                eventRate_grad(1) = Tb * nBackgrPerDay;
            end
        end
        % by A
        eventRate_grad(2)   = Tb * triggRate/A;
        % by alpha
        eventRate_grad(3)   = Tb * sum( triggRates_ij .* (mag(isContributing)-Mc) );
        % by c
        tmp_c               = -p ./ fTempK_inn_ij; % * ci(isContributing)
        eventRate_grad(4)   = Tb * sum( triggRates_ij .* tmp_c );
        % by p
        tmp_p               = -log(fTempK_inn_ij);
        eventRate_grad(5)   = Tb * sum( triggRates_ij .* tmp_p );
        % by D, gamma, q
        if ismember(strModus, {'time-space rates', 'decluster'})
            helpTerm_D          = (q * (1-fSpatK_inn_ij.^(-1)) - 1 - norm_D(isContributing)); % .* isWiderThanMinWidth(isContributing);
            helpTerm_q          = (1/(q-1) - log(fSpatK_inn_ij) - norm_q(isContributing));
            eventRate_grad(6)   = sum( triggRates_ij .* 1/D .* helpTerm_D );
            eventRate_grad(7)   = sum( triggRates_ij .* (mag(isContributing)-Mc) .* helpTerm_D );
            eventRate_grad(8)   = sum( triggRates_ij .* helpTerm_q );
        elseif strcmp( strModus, 'time rates')
            eventRate_grad(6)   = Tb * sum( factor_dD(isContributing) .* fTempK_inn_ij_toMinusP );
            eventRate_grad(7)   = Tb * sum( factor_dGamma(isContributing) .* fTempK_inn_ij_toMinusP );
            eventRate_grad(8)   = Tb * sum( factor_dQ(isContributing) .* fTempK_inn_ij_toMinusP );
        end
        % by Tb
        if strcmp( strModus, 'time rates')
            eventRate_grad(9)   = eventRate/Tb;
        end
        
        % Overwrite gradients for fixed parameters by 0
        eventRate_grad(fixedParamETAS) = 0;
    
    end
                                                     
end