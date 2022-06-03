function tAnalysis = create_analysisTable4catalogs( nRows )

    tAnalysis                   = table;
    tAnalysis.catalogID         = zeros(nRows,1);
    tAnalysis.nEvents           = zeros(nRows,1);
    tAnalysis.nTargetEvents     = zeros(nRows,1);
    tAnalysis.percTriggered     = zeros(nRows,1);
    tAnalysis.largestMagn       = zeros(nRows,1); 
    tAnalysis.secondMagn        = zeros(nRows,1); 
    tAnalysis.nLarger5          = zeros(nRows,1); 
    tAnalysis.tempVolatility    = zeros(nRows,1); 
    tAnalysis.spatVolatility    = cell(nRows,1);
    tAnalysis.spatDistr         = cell(nRows,1);
    tAnalysis.doubletSampleSize = cell(nRows,1);
    tAnalysis.doubletPercentage = cell(nRows,1);
    tAnalysis.doubletType       = cell(nRows,1);
    tAnalysis.bathLawMagnDiff   = cell(nRows,1);
    tAnalysis.bathLawType       = cell(nRows,1);
    
end