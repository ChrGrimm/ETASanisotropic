function mu_backgrPerDay = convert_mu2backgrPerDay( Catalog, Inputs, SpatData, mu )

    if strcmp(Inputs.ModelSettings.spaceModel, 'none')
        mu_backgrPerDay     = mu;
        
    else
        isTargetT       = Catalog.flag >= 0;
        tWindow         = Inputs.TargetSettings.tWindow_days;
        length_tWindow  = sum( tWindow(:,2)-tWindow(:,1) );
        mu_backgrPerDay = mu * sum( Catalog.backgrProb(isTargetT) .* SpatData.backgrIntegral ) / length_tWindow;
        
    end
    
end