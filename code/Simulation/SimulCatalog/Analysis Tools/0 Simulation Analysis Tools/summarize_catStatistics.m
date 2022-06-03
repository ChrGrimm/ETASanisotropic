function CatStatistics = summarize_catStatistics( CatStatistics, iRun )

    load('Simulation_results1.mat','tAnalysis_synthCat')
    load('Estimation_results.mat','ModelAnalysis') 
    
    CatStatistics.nTargetEvents(:,iRun)   = tAnalysis_synthCat.nTargetEvents;
    CatStatistics.largestMagn(:,iRun)     = tAnalysis_synthCat.largestMagn;
    CatStatistics.tempVolatility(:,iRun)  = tAnalysis_synthCat.tempVolatility;
    CatStatistics.spatVolatility(:,iRun)  = tAnalysis_synthCat.spatVolatility(1:1000);
    CatStatistics.integrRateTempVol(iRun) = ModelAnalysis.tempVolatility;

end