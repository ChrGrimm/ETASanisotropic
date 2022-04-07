function estimate_etasModel(pathIniFile, pathResultsFile, makeLog, saveDiagnosticPlots )
% This is the main function of the ETAS estimation process. Call this function, if you want to
% estimate an ETAS model for given input settings defined in pathIniFile (.mat file).
% The function stores the estimation results in a .mat file in pathResultsFile.
% 
% INPUTS:
% - pathIniFile:         string     path to .mat ini file with ETAS estimation setting structures
% - pathResultsFile:     string     path to .mat results file where estimation results are saved
%                                   (if path does not yet exist, directory is created)
% - makeLog:             boolean    If true, a logfile is saved in the results folder (folder path of pathResultsFile)
% - saveDiagnosticPlots: boolean    If true, diagnostic plots are created and saved in the results folder (folder path of pathResultsFile)
%
    
    %% Check pathResultsFile and, if needed, create directory
    [pathResultsFolder, resultsFilename] = check_pathResultsFile( pathResultsFile );
    
    %% Start logfile
    tStart = tic;
    if makeLog
        startLogFile(pathResultsFolder, resultsFilename)
    end
    
    disp( ['... Run ETAS Estmation of ', pathResultsFolder, '...'] )
    
    %% Preprocess inputs: Compile catalog, precompute spatial and temporal data
    % Load ini file
    EstimSettings = load(pathIniFile);
    % Preprocess inputs
    [ Catalog, ...
      ExtendedSettings, ...
      ModelFuncs, ...
      SpatData, ...
      TempData ] = preprocess_estimationInputs( EstimSettings, makeLog );
    
    %% ETAS parameter estimation: Parameters are iteratively optimized by Log-Likelihood estimation
    [ Catalog, ...
       IterationLog, ...
       paramETAS, ...
       SpatData, ...
       triggProbMatrix ] = fit_etasModel( Catalog, ExtendedSettings, ModelFuncs, SpatData, TempData, makeLog );
    
    %% Postprocess results, compute model summaries
    [ ModelSummary, ...
      ModelFuncs, ...
      paramETAS_descrip, ...
      TriggerRelationsETAS, ...
      ClustersETAS ] = postprocess_estimationResults( Catalog, ...
                                                        ExtendedSettings, ...
                                                        ModelFuncs, ...
                                                        SpatData, ...
                                                        paramETAS, ...
                                                        triggProbMatrix, ...
                                                        IterationLog );
    
    IterationLog.elapsedTime = toc(tStart)/60;
    
    %% Print some results                                                   
    % Performance time   
    disp(['ETAS model estimation took ', num2str(IterationLog.elapsedTime), ' minutes.']) 
    % b value estimate
    disp([ 'Model estimates b value = ', num2str(paramETAS(10)/log(10)) ])
    % Integrated ETAS rate over entire time-space window.
    % Should closely approximate the number of target events (flag>0, without duplicates) in the original catalog.
    disp(['Models estimates ', num2str(IterationLog.integratedRate(end)), ' target events (original catalog: ', num2str(sum(Catalog.flag>0 & ~Catalog.isDupl)), ' events)'])
    
    %% Store results in .mat file in model results folder
    diary off
    save(pathResultsFile, 'paramETAS', 'paramETAS_descrip', 'IterationLog', 'ModelSummary', 'EstimSettings', 'ModelFuncs', 'ClustersETAS'); % 'TriggerRelationsETAS', 
    
    %% Store some diagnostic plots
    % Plot observed vs. integrated rate
    figure
    plot_etasRateVsObservedEvents( ExtendedSettings, Catalog, paramETAS, SpatData, TempData, ModelFuncs, ModelSummary.mu_backgrPerDay, pathResultsFolder );
    savefig(fullfile(pathResultsFolder, ['Plot_etasVsObservedRates_', resultsFilename, '.fig']))
    
    if saveDiagnosticPlots
        create_estimDiagnosisPlots( pathResultsFile, pathResultsFolder, resultsFilename, ExtendedSettings )
    end

end