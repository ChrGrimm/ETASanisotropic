function Catalog = estimate_nodalPlaneSolutions( Catalog, ...
                                                origCatalog, ...
                                                Inputs, ...
                                                ModelFuncs, ...
                                                idxAniso )
   
    %% Estimation by USGS slab model 2.0    
    if strcmp(Inputs.SpaceSettings.sampleStrikes, 'USGS slab model 2.0')
        % Load USGS strike data and compute scattered interpolant
        load('strike_usgsSlab2.mat','Strike_usgsSlab2')
        F_strike = scatteredInterpolant(Strike_usgsSlab2(:,1), Strike_usgsSlab2(:,2), ...
                                        Strike_usgsSlab2(:,3), 'nearest', 'nearest');
        
        % Interpolate strike
        Catalog.strike(idxAniso)  = F_strike(Catalog.lon(idxAniso), Catalog.lat(idxAniso));
        Catalog.epiPos(idxAniso)  = 0.5*ones(length(idxAniso),1);
    
    %% Estimation via forward trigger rate maximization
    else
    % strcmp(SpatKernel.methodStrikeEstim, 'max trigger rate')
        
        if isnumeric(Inputs.SpaceSettings.sampleStrikes)
    %         strikeSample        = [ origCatalog.str1(idxAniso), origCatalog.str2(idxAniso), ones(length(idxAniso),1) * Inputs.SpaceSettings.sampleStrikes(:)' ];
            strikeSample        = ones(length(idxAniso),1) * Inputs.SpaceSettings.sampleStrikes;
            epiPosSample        = Inputs.SpaceSettings.sampleEpicenter; 
            
        elseif ismember( SpatKernel.sampleStrikes, {'gcmt', 'iscgem', 'self'} )
            % - match nodal plane solutions from external catalogs (gcmt or iscgem)
            % Enrich (strike-dip-rake) sets from external dataset
            origCatalog = match_nodalPlaneSolutions( origCatalog, ...
                                                     Inputs.TargetSettings.polygonXY, ...
                                                     Inputs.SpatialSettings.sampleStrikes );
            % Set rupture length = 0 (numerically) for events that do not have nodal plane provided                                            
            isSetsGiven         = ~(origCatalog.str1==0 & origCatalog.str2==0);
            % Define strike sample based on nodal plane sets as enriched before
            strikeSample    = [ origCatalog.str1( isSetsGiven ), origCatalog.str2( isSetsGiven ) ];
            
        end
        
        % - choose the set that is providing higher trigger rate under initial parameter guesses
        [ Catalog.strike(idxAniso), ...
          Catalog.epiPos(idxAniso) ]  = estimate_byMaxTriggerDensity( strikeSample, ...
                                                                        epiPosSample, ...
                                                                        Catalog, ...
                                                                        Inputs, ...
                                                                        ModelFuncs, ...
                                                                        idxAniso );
                                                                    
        % Events that are better modelled isotropically
        overwriteByIso                      = Catalog.strike==-1;
        Catalog.typeKernel(overwriteByIso)  = {'iso'};
        Catalog.rupExtent( overwriteByIso ) = 0;
        Catalog.strike( overwriteByIso )    = NaN;
        Catalog.epiPos( overwriteByIso )    = NaN;
        Catalog.spatRestr(overwriteByIso)   = ModelFuncs.fRestrSpace_degrees( Catalog.mag(overwriteByIso), 'iso' );
    
    end
    
    %% Checks     
    % analyze_strikeDifferences( origCatalog, Catalog.strike );
    % plot_rotatedAftershockClouds
       
end
