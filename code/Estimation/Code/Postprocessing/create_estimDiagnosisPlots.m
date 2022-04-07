function create_estimDiagnosisPlots( pathResultsFile, pathResultsFolder, resultsFilename, EstimSettings )

    % Estimated aftershock productivity (direct aftershocks only)
    figure 
    plot_aftershockProductivity( pathResultsFile, EstimSettings.TimeSettings.restr_days, 'b', '-', resultsFilename, true )
    savefig(fullfile(pathResultsFolder, ['Plot_AftershProd_', resultsFilename, '.fig']))

    % Estimated cluster sizes (direct and secondary aftershocks)
    figure
    plot_clusterSize( pathResultsFile, EstimSettings.TimeSettings.restr_days, 'b', '-', resultsFilename, true )
    savefig(fullfile(pathResultsFolder, ['Plot_ClusterSize_', resultsFilename, '.fig']))

    % Omori Law
    figure
    plot_omoriLaw( pathResultsFile, 7, 'b', '-', resultsFilename )
    savefig(fullfile(pathResultsFolder, ['Plot_OmoriLaw_', resultsFilename, '.fig']))

    % 3D pdf spatial kernel
    figure
    plot_pdfSpatialKernel( pathResultsFile, [5,6,7] )
    savefig(fullfile(pathResultsFolder, ['Plot_pdfSpatialKernel_', resultsFilename, '.fig']))

    % cdf spatial kernel
    figure
    plot_cdfSpatialKernel( pathResultsFile, [5,6,7], 'b', '-', resultsFilename )
    savefig(fullfile(pathResultsFolder, ['Plot_cdfSpatialKernel_', resultsFilename, '.fig']))

end