function [id, depth, strike, epiPos] = match_aftershocks2met( xSampled, ...
                                                              ySampled, ...
                                                              magSampled, ...
                                                              MET )
   
    %% Preallocate output
    nAftershocks    = length(xSampled);
    id              = zeros(nAftershocks, 1);
    depth           = zeros(nAftershocks, 1);
    strike          = zeros(nAftershocks, 1);
    epiPos          = zeros(nAftershocks, 1);
    
    %% Loop over magnitudes
    eps                 = 10^-5;
    for Mw = unique(magSampled)'
        idxAftersh          = find( abs(magSampled-Mw) < eps );
        idxMET              = find( abs(MET.mag-Mw) < eps );
        % Find closest MET location to sampled x and y
        idxNearest = knnsearch( [MET.x(idxMET),MET.y(idxMET)], [xSampled(idxAftersh),ySampled(idxAftersh)] );

        idxMatch            = idxMET(idxNearest);                     
        id(idxAftersh)      = MET.eventID(idxMatch);
        depth(idxAftersh)   = MET.depth(idxMatch);
        strike(idxAftersh)  = MET.strike(idxMatch);
        epiPos(idxAftersh)  = MET.epiPos(idxMatch);
    end
    
end