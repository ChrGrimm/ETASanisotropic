%% Manual Inputs
% Model path (from folder Results)
Models = {'CAL-Sequences\2019_RidgecrestMw64_ETASI_220111_aniso'};
      
% Model names to appear in figure legends      
LegendNames = {'test'};

% Choose type of analysis (choose from cases below) or 'none' to only
% obtain table with ETAS model result statistics
typeAnalysis = 'productivity';

%% Code
figure
CatStatistics   = struct;

%% Create plots and statistics
nModels = length(Models);

for iRun = 1:nModels

    % Add directory to path temporarily
    manage_etasModelRunPaths( Models{iRun}, 'add' )

    % Table with measures 
    tModelComparison(iRun,:) = create_tableOfEtasResults( Models{iRun}, LegendNames{iRun} );

    switch typeAnalysis
        case 'nAftershocks'
            cdfPlot_numberOfAftershocks( LegendNames{iRun}, iRun, iRun==nModels )

        case 'largestAftersh'
            pdfPlot_magnOfStrongestAftersh( LegendNames{iRun}, iRun, iRun==nModels )

        case 'spatDistr'
            scatterPlot_aftershLocations

        case 'productivity'
            plot_aftershProductivity( LegendNames{iRun}, iRun, typeAnalysis )

        case 'cluster size'
            plot_clusterSize( LegendNames, iRun )

        case 'spatial cdf'
            plot_spatialCDF( LegendNames, iRun )

        case 'spatial pdf 3D'
            plot_3dSpatialPDF( LegendNames, iRun )

        case 'bath law'
            plot_bathLaw( LegendNames, iRun )

        case 'doublets'
            plot_doublets( LegendNames, iRun )

        case 'type of aftershock'
            plot_typeOfAftershock( LegendNames, iRun )

        case 'catalog statistics'
            CatStatistics = summarize_catStatistics( CatStatistics, iRun );

        case 'monthly rates'
            load('Estimation_results.mat','Catalog','TargetWindow','Method','paramETAS','ModelAnalysis') 
            dailyBackgrRate = ModelAnalysis.backgrRatePerDay;
            figure
            compute_aggregEventRate( Catalog, TargetWindow, Method, sqrt(paramETAS), dailyBackgrRate, true );
%                 title(LegendNames{iRun})

    end

    rmpath(genpath(modelRunDirectory))

end

switch typeAnalysis
    case 'doublets'
%             plot_doubletsOrigCatalogs( cellStrRuns{iRun} )

    case 'catalog statistics'
        plot_catalogStatistics( CatStatistics, LegendNames )

    case 'bath law'
        % Ideal Bath Law line
        plot( [5.5; 9.0] , [1.2; 1.2], 'k:', 'LineWidth', 1.5, 'DisplayName', 'Bath`s Law')

    case 'latexResults'
        create_latexResultsTable( tModelComparison );

end
