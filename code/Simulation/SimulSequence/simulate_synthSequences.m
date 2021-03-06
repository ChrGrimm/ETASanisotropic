function [ tAnalysis, ...
           Analysis_detail ] = simulate_synthSequences( Settings4SynthSeq, ...
                                                           DoubletCriteria, ...
                                                           Catalog, ...
                                                           TargetWindow, ...
                                                           SpatKernel, ...
                                                           Tectonics, ...
                                                           GeneralSettings, ...
                                                           OmoriLaw, ...
                                                           paramETAS, ...
                                                           beta, ...
                                                           ModelAnalysis )
    
    tStart = tic;
                                                       
    %% Simulation & Analysis
    rng(1)
    iSetSimul = 0;
    
    %% Prepare simulation
    % Initial catalog
    [ Settings4SynthSeq, ...
      magnDistr, ...
      productivity, ...
      tAnalysis, ...
      Analysis_detail ] = prepare_synthSequenceSimulation( Settings4SynthSeq, ...
                                                            ModelAnalysis, ...
                                                            Catalog, ...
                                                            TargetWindow, ...
                                                            SpatKernel, ...
                                                            Tectonics, ...
                                                            OmoriLaw, ...
                                                            paramETAS, ...
                                                            beta );
   
    % Loop over initial trigger magnitudes
    for iTriggerMagn = Settings4SynthSeq.triggerMagn
        
        iSetSimul = iSetSimul + 1;
        
        %% Prepare trigger event
        rupL                    = estimate_rupSize( iTriggerMagn, Tectonics.type, Tectonics.faultingStyle, GeneralSettings.spaceUnit );
        [triggerX, triggerY]    = centroid( polyshape( TargetWindow.polygonXY(:,1), TargetWindow.polygonXY(:,2) ) );
        triggerEvent            = array2table( [ 1, 0, triggerX, triggerY, iTriggerMagn, -1, rupL, 0, -1, 1 ], ...
                                               'VariableNames',{'id','t','x','y','mag','flag','rupL','strike','triggerID','clusterID'} );
                                                         
        %% Loop over realizations   
        SynthSequences = cell(Settings4SynthSeq.nRealizations,1);
        
        for iRealiz = 1:Settings4SynthSeq.nRealizations

            %% Produce synthetic catalogs
            SynthSequences{iRealiz} = simulate_synthEventSets( triggerEvent, ...
                                                                magnDistr, ...
                                                                paramETAS, ...
                                                                productivity, ...
                                                                Settings4SynthSeq, ...
                                                                SpatKernel, ...
                                                                Tectonics, ...
                                                                GeneralSettings, ...
                                                                TargetWindow, ...
                                                                ModelAnalysis.f_inn );
                                                            
        end
        
        %% Analyze catalogs   
        tAnalysis.triggerMagn{iSetSimul} = ['Mw', num2str(10*iTriggerMagn)];
        [ Analysis_detail(iSetSimul), ...
          tAnalysis(iSetSimul,2:end) ] = analyze_sequences( SynthSequences, ...
                                                             DoubletCriteria, ...
                                                             ['Mw', num2str(10*iTriggerMagn)], ...
                                                             GeneralSettings.spaceUnit, ...
                                                             TargetWindow.M_c );
                                                     
    end
    
    tEnd = toc(tStart)/60;
    disp(['Elapsed time for sequence simulation is ', num2str(tEnd), ' minutes.'])
        
end