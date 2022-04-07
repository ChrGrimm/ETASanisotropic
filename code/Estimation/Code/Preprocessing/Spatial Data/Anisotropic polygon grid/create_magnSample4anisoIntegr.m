function [ sample4anisoIntegr, ...
           anisoPeaks ] = create_magnSample4anisoIntegr( magnitudes, ...
                                                         ModelFuncs )

    uniqueMagn              = unique(magnitudes);
    sample4anisoIntegr      = zeros(length(uniqueMagn), 3);
    sample4anisoIntegr(:,1) = uniqueMagn; % unique(Catalog.mag(Catalog.rupExtent>0)); % (0:0.1:max(Catalog.mag))';
    sample4anisoIntegr(:,2) = ModelFuncs.fRupLength( uniqueMagn ); % estimate_rupSize(uniqueMagn, SpaceSettings, Tectonics.faultingStyle, spaceUnit ); 
    sample4anisoIntegr(:,3) = ModelFuncs.fRestrSpace_degrees( uniqueMagn, 'aniso'); % max(k * sample4anisoIntegr(:,2), c);
    
    anisoPeaks(:,1)         = sample4anisoIntegr(:,3);
    anisoPeaks(:,2)         = min(0.1, sample4anisoIntegr(:,2));

end