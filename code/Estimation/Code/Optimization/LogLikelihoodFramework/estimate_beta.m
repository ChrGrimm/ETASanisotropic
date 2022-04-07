function beta_llh = estimate_beta( Catalog, Mc )

    %% Loglikelihood 
    isTarget = Catalog.flag==1;
    % Estimate according to Jalilian, 2019
    beta_llh = sum(isTarget) / sum(Catalog.mag(isTarget)-Mc); 
    
%     %% Root mean square
%     % Count number of events exceeding each magnitude level
%     mag             = ( Mc : 0.1 : round(Mc+max(Catalog.mag(isTarget)),1) )';
%     nEvents         = hist(round(Mc+Catalog.mag(isTarget),1), mag)';
%     % Determine annual occurrence rates
%     rateEvents      = nEvents/round(diff(timeWindow_days)/365);
%     rateEventsCum   = cumsum(rateEvents, 'reverse');
%     % Compute linear regression
%     A = [ones(length(mag),1), mag];    % Vandermonde matrix
%     f = log10(rateEventsCum);           % data vector (right hand side)
%     c = A\f;                            % coefficients of overdetermined system
%     beta_rms = -c(2)*log(10);
    
end