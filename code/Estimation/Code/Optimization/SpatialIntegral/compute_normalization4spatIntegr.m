function [ normaliz, ...
           normalizSqu, ...
           normaliz_D, ...
           normaliz_q ] = compute_normalization4spatIntegr( Catalog, ...
                                                            isIso, ...
                                                            fSpatK_inn, ...
                                                            q )

    % Extract columns
    R                   = Catalog.spatRestr(isIso);
    mag                 = Catalog.mag(isIso);
    rupExtent           = Catalog.rupExtent(isIso);

    % Evaluate inner spatial function until spatial extent
    f_inn_R             = fSpatK_inn(R, R.^2, mag, rupExtent);

    % Compute normalization as integral until spatial extent
    f_inn_R_to1minusQ   = f_inn_R.^(1-q);
    normaliz            = 1 - f_inn_R_to1minusQ;
    normalizSqu         = normaliz.^2;

    % Compute normalizations of gradient integrals
    normaliz_D          = f_inn_R.^(-q) .* (f_inn_R - 1);
    normaliz_q          = log(f_inn_R) .* f_inn_R_to1minusQ;
        
    end