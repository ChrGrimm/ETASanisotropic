Models = {PUT IN LIST OF MODEL RESULT FOLDERS WHERE YOU WANT TO MODIFY INPUTS};

%% Loop over models
for iModel = 1:length(Models)
    
    cd(Models{iModel})
    clear all
    
    %% Load inputs .mat file, modify it, save it
    
end