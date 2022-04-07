function [pathResultsFolder, resultsFilename] = check_pathResultsFile( pathResultsFile )

    [pathResultsFolder, resultsFilename, fileType] = fileparts(pathResultsFile);
    if ~strcmp(fileType, '.mat')
        error('pathResultsFile must have ''.mat'' ending.')
    end
    
    if ~exist(pathResultsFolder, 'dir')
        mkdir(pathResultsFolder)
    elseif exist(pathResultsFile, 'file')
        error(['File named ', resultsFilename, ' already exists in directory ', pathResultsFolder, '.'])
    end

end