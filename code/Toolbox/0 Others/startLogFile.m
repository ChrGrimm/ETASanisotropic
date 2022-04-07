function startLogFile(pathResultsFolder, resultsFilename)

    logFilename = ['Logfile_', resultsFilename, '.log'];
    pathLogfile = fullfile(pathResultsFolder, logFilename);
    if exist(pathLogfile, 'file')
        diary off
        error(['File named ', logFilename, ' already exists in directory ', pathResultsFolder, '. Consider removing or renaming it.'])
    else
        diary(pathLogfile)
    end    

end