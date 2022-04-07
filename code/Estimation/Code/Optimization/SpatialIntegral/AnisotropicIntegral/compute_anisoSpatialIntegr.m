function [ F, ...
           dF_D, ...
           dF_q, ...
           anisoPeaks ] = compute_anisoSpatialIntegr( inclGradients, ...
                                                      Catalog, ...
                                                      isAniso, ...
                                                      spatParamETAS, ...
                                                      fSpatK_inn, ...
                                                      fSpatK_factor, ...
                                                      SpatData )
       
    %% Anisotropic spatial integral
    % Extract (isotropic) event-to-polygon distances
    r_iso               = SpatData.Iso_r( isAniso );
    r0_iso              = SpatData.Iso_r0( isAniso );
    segmFactor_iso      = SpatData.Iso_segmFactor( isAniso );
    isStartP_iso        = SpatData.Iso_isStartP( isAniso );
    isEndP_iso          = SpatData.Iso_isEndP( isAniso );
    % Extract (anisotropic) rupture-to-polygon distances, segment weights and markers
    r_aniso             = SpatData.Aniso_r( isAniso );
    r0_aniso            = SpatData.Aniso_r0( isAniso );
    segmFactor_aniso    = SpatData.Aniso_segmFactor( isAniso );
    isStartP_aniso      = SpatData.Aniso_isStartP( isAniso );
    isEndP_aniso        = SpatData.Aniso_isEndP( isAniso );
    
    mag                 = Catalog.mag( isAniso );

    % Compute interpolation sets for isotropic and anisotropic integral parts
    [ sample_integrAniso, ...
      sample_integrIso, ...
      sample_gradAniso_dD, ...
      sample_gradIso_dD, ...
      sample_gradAniso_dQ, ...
      sample_gradIso_dQ, ...
      distGrid_integr, ...
      distGrid_grad, ...
      anisoPeaks ] = compute_anisoIntegrSamples( inclGradients, ...
                                                 spatParamETAS, ...
                                                 SpatData.magnSample4anisoIntegr, ...
                                                 SpatData.anisoPeaks, ...
                                                 fSpatK_inn, ...
                                                 fSpatK_factor );

    nEvents = sum(isAniso);
    % Initialize output variables
    [ F, dF_D, dF_q ]   = initialize_spatIntegrOutputVectors( Catalog(isAniso,:), 'default' );
        
    % Loop over all events
    for iEv = 1:nEvents

        % Extract relevant samples
        magnIndex       = find(SpatData.magnSample4anisoIntegr(:,1)==mag(iEv)); % 10*round(magn_iEv,1)+1;

        % Interpolate anisotropic integral and gradiants
        [ integralAniso, ...
          gradiantAniso_D, ...
          gradiantAniso_q ]= interpolate_integrals( sample_integrAniso(:,magnIndex), ...
                                                    sample_gradAniso_dD(:,magnIndex), ...
                                                    sample_gradAniso_dQ(:,magnIndex), ...
                                                    distGrid_integr(:,magnIndex), ...
                                                    distGrid_grad(:,magnIndex), ...
                                                    inclGradients, ...
                                                    r_aniso{iEv}, ...
                                                    r0_aniso{iEv}, ...
                                                    segmFactor_aniso{iEv}, ...
                                                    isStartP_aniso{iEv}, ...
                                                    isEndP_aniso{iEv} );

        % Interpolate isotropic integral and gradiants
        [ integralIso, ...
          gradiantIso_D, ...
          gradiantIso_q ]= interpolate_integrals( sample_integrIso(:,magnIndex), ...
                                                    sample_gradIso_dD(:,magnIndex), ...
                                                    sample_gradIso_dQ(:,magnIndex), ...
                                                    distGrid_integr(:,magnIndex), ...
                                                    distGrid_grad(:,magnIndex), ...
                                                    inclGradients, ...
                                                    r_iso{iEv}, ...
                                                    r0_iso{iEv}, ...
                                                    segmFactor_iso{iEv}, ...
                                                    isStartP_iso{iEv}, ...
                                                    isEndP_iso{iEv} );                                        

        % Sum up anisotropic and isotropic contributions
        F(iEv)    = integralAniso + integralIso;
        dF_D(iEv) = gradiantAniso_D + gradiantIso_D;
        dF_q(iEv) = gradiantAniso_q + gradiantIso_q;

    end
    
end
                  