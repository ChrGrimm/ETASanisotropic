function SpatData = initialize_tables4precomputations( SpatData, nEvents )
       
    SpatData.Iso_r              = cell(nEvents, 1);
    SpatData.Iso_r0             = cell(nEvents, 1);
    SpatData.Iso_r_squ          = cell(nEvents, 1);
    SpatData.Iso_r0_squ         = cell(nEvents, 1);
    SpatData.Iso_segmFactor     = cell(nEvents, 1);
    SpatData.Iso_isStartP       = cell(nEvents, 1);
    SpatData.Iso_isEndP         = cell(nEvents, 1);

    SpatData.Aniso_r                = cell(nEvents, 1);
    SpatData.Aniso_r0               = cell(nEvents, 1);
    SpatData.Aniso_r_squ            = cell(nEvents, 1);
    SpatData.Aniso_r0_squ           = cell(nEvents, 1);
    SpatData.Aniso_segmFactor       = cell(nEvents, 1);
    SpatData.Aniso_isStartP         = cell(nEvents, 1);
    SpatData.Aniso_isEndP           = cell(nEvents, 1);
    
    SpatData.magnSample4anisoIntegr = [];
    SpatData.anisoPeaks             = [];
    
end