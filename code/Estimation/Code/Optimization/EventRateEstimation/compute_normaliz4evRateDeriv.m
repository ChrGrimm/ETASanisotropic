function [ norm_D, norm_q ] = compute_normaliz4evRateDeriv( fSpatK_inn, R, mag, rupExtent, q )
                                                                      
        spatK_inn_R             = fSpatK_inn(R, R.^2, mag, rupExtent);
        spatK_inn_R_to1minusQ   = spatK_inn_R.^(1-q);
        normalization           = 1 - spatK_inn_R_to1minusQ;
        norm_D                  = (1-q) * (spatK_inn_R-1) ./ (normalization .* spatK_inn_R.^q);
        norm_q                  = log(spatK_inn_R) .* spatK_inn_R_to1minusQ ./ normalization;  
        
end