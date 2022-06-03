function plot_bathLaw( legendNames, iRun )

    %% Version
    strVersion  = 'shared plot';

    %% Load data
    load('Estimation_results.mat','TargetWindow')
    load('Simulation_results.mat','tAnalysis_synthSeq','tAnalysis_synthCat')
%     load('Simulation_resultsMainshocksOnly.mat','tAnalysis_synthSeq','tAnalysis_synthCat')
    % Bath law data for original catalog
%     load('Results_histCatalogs.mat')
    load('Results_histCats_mainshockOnly.mat')
    if strcmp(TargetWindow.region, 'CAL')
        bathLaw_historic = tAnalysis_historicCat.bathLawMagnDiff{2};      
    elseif strcmp(TargetWindow.region, 'JPN')
        bathLaw_historic = tAnalysis_historicCat.bathLawMagnDiff{1};
    end
    
    magnitudes = (5.5:0.1:max(tAnalysis_synthSeq.mainMw))';
    
%     strLegend   = legendNames{ iRun };
    
    %% Preprocess bath law results to enable plotting
    [ bathLaw_seq1, ...
      bathLaw_seq2, ...
      bathLaw_cat ] = preprocess_bathLawResults( tAnalysis_synthSeq, ...
                                                tAnalysis_synthCat );
    
    %% Single Plot
    if strcmp(strVersion, 'single plot')
        singleplot_bathLaw( magnitudes, bathLaw_seq1, bathLaw_cat, bathLaw_historic, legendNames{iRun} )
    
    %% Shared plot
    elseif strcmp(strVersion, 'shared plot')
        sharedplot_bathLaw( magnitudes, bathLaw_seq1, bathLaw_cat, bathLaw_historic, legendNames{iRun}, TargetWindow.region, iRun )
                
    else
        error('Unknown plotting version')
    end
    
    t = title('(a)', 'FontSize', 10);
    t.Units = 'Normalize'; 
    t.Position(1) = 0;
    t.HorizontalAlignment = 'left';
    
end