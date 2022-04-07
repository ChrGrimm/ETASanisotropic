function [ F, ...
            dF_D, ...
            dF_gamma, ...
            dF_q ] = initialize_spatIntegrOutputVectors( Catalog, spaceModel )

    if strcmp(spaceModel, 'none')
        % Spatial integrals are set to 1, gradiants to 0
        F           = abs(Catalog.flag) >= 0.5;
        dF_D        = 0;
        dF_gamma    = 0;
        dF_q        = 0;
        
    else
        nEvents     = size(Catalog,1);
        % Initialize output variables
        F           = zeros(nEvents,1);
        dF_D        = zeros(nEvents,1);
        dF_gamma    = zeros(nEvents,1);
        dF_q        = zeros(nEvents,1);
        
    end
        
end