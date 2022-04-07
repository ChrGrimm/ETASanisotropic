function [ exactEstim_N0, ...
            exactEstim_dN0 ] = test_N0_integral( Catalog, ...
                                                   timeGrid, ...
                                                   M_c, ...
                                                   nBackgrPerDay, ...
                                                   sqrtParamETAS, ...
                                                   F, ...
                                                   productivity, ...
                                                   strMode )

    computeGradients    = strcmp(strMode, 'with gradient');
    
    % Extract ETAS parameters (they are in square root format  for numerical reasons)
    [mu, A, ~, c, p, ~, ~, ~, Tb] = extract_paramETASvalues( sqrtParamETAS, 'sqrt' );
    
    % Extract event-specific information in columns of table Catalog
    t           = Catalog.t;
    mag         = Catalog.mag;
    evWeight    = Catalog.evWeight;
    
    factor          = productivity .* evWeight .* F;
%     factor_dD       = productivity .* evWeight .* dF_D;
%     factor_dGamma   = productivity .* evWeight .* dF_gamma;
%     factor_dQ       = productivity .* evWeight .* dF_q;
    
    idxTargetTime   = find(Catalog.flag>=0);
    exactEstim_N0   = zeros(1, length(timeGrid));
    exactEstim_dN0  = zeros(10, length(timeGrid));
    
    for jCounter = 1:length(timeGrid)
        
        if jCounter==31069
            test = 1;
            exact_isContribEvent = isContribEv;
            exact_idxContribEvent = find(isContribEv);
            exact_tContr = tContr;
            exact_factorContr = factorContr;
            exact_trigij = funcN0_trig( timeGrid(jCounter) );
            exact_N0 = funcN0( timeGrid(jCounter) );
        end
                
        isContribEv = t <= timeGrid(jCounter); % & t >= tStart-0.1;
        tContr      = t(isContribEv);
        factorContr = factor(isContribEv);
%         if ~isempty(tContr)
%             if length(tContr)==1
%                 tContr      = [tContr; tContr];
%                 factorContr = [factorContr; 0];
%             end
            funcN0_trig = @(tt) factorContr .* (tt-tContr+c).^(-p);
            funcN0      = @(tt) Tb * ( mu * nBackgrPerDay + sum( funcN0_trig(tt) ) );
            funcLL2      = @(tt) (1 - exp(- funcN0(tt) )) / Tb;
            funcLL2_mu   = @(tt) exp(- funcN0(tt) ) * nBackgrPerDay * 2 * sqrtParamETAS(1);
            funcLL2_A    = @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt)/A ) * 2 * sqrtParamETAS(2);
            funcLL2_alpha= @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt) .* (mag(isContribEv)-M_c) ) * 2 * sqrtParamETAS(3);
            funcLL2_c    = @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt) .* (-p./(tt-tContr+c)) ) * 2 * sqrtParamETAS(4);
            funcLL2_p    = @(tt) exp(- funcN0(tt) ) .* sum( funcN0_trig(tt) .* (-log(tt-tContr+c)) ) * 2 * sqrtParamETAS(5);
            funcLL2_Tb   = @(tt) (-1/Tb^2) + exp(- funcN0(tt) ) / Tb .* ( mu * nBackgrPerDay + sum( funcN0_trig(tt) ) + 1) * 2 * sqrtParamETAS(9);

            exactEstim_N0(jCounter)     = funcLL2( timeGrid(jCounter) );
            exactEstim_dN0(1,jCounter)  = funcLL2_mu( timeGrid(jCounter) );
            exactEstim_dN0(2,jCounter)  = funcLL2_A( timeGrid(jCounter) );
            exactEstim_dN0(3,jCounter)  = funcLL2_alpha( timeGrid(jCounter) );
            exactEstim_dN0(4,jCounter)  = funcLL2_c( timeGrid(jCounter) );
            exactEstim_dN0(5,jCounter)  = funcLL2_p( timeGrid(jCounter) );
            exactEstim_dN0(9,jCounter)  = funcLL2_Tb( timeGrid(jCounter) );
            
%         end
    end
