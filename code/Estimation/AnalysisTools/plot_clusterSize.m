function plot_clusterSize( pathResultsFile, nDays, color, lineStyle, legendEntry, isLogScale )
% This function plots the expected cluster size against triggering magnitude, given an ETAS model 
% results .mat file. This plot is similar to aftershock productivity, but accounts for secondary 
% triggerimg, also.
%
% INPUTS:
% - pathResultsFile:        string      full path to .mat file containing ETAS estimation results
% - nDays:                  scalar      duration (in days) of aftershock triggering (e.g. 10, 100, 365)
% - color:                              line color code (Matlab notation)
% - lineStyle:                          line style code (Matlab notation)
% - legendEntry:            string      name of model to appear in legend
% - isLogScale:             boolean     If true, y axis is in log scale (recommended)

    %% Load or set ETAS parameters
    if ~isempty(pathResultsFile)
        % Load parameter estimates from submitted path to results file
        load(pathResultsFile, 'paramETAS', 'EstimSettings', 'ModelSummary')
        [~, A, alpha, c, p] = extract_paramETASvalues( paramETAS, 'default' );
        Mc                  = EstimSettings.TargetSettings.Mc;
        Mmax                = EstimSettings.TargetSettings.Mmax;
        branchingRatio      = ModelSummary.branchRatio(1);
    else
        % Manually set parameter estimates
        A               = 0.01;
        alpha           = 2;
        c               = 0.01;
        p               = 0.9;
        Mc              = 4;
        Mmax            = 9;
        branchingRatio  = 0.5;
    end
    
    if branchingRatio > 1
        warning('Cluster size is expected to be infinity, due to branching ratio > 1.')
        
    else
    
        %% Compute cluster size function
        % As the Omori law is not normalized (does not integrate to 1), scale aftershock productivity
        OmoriIntegral       = 1/(1-p) * ( (nDays+c).^(1-p) - c.^(1-p) );
        funcClusterSize     = @(mag, Mc) OmoriIntegral * A / (1-branchingRatio) * exp(alpha*(mag-Mc));
        % Magnitude vector
        magnitudes = Mc:0.1:Mmax;

        %% Plotting
        hold on
        plot( magnitudes, funcClusterSize(magnitudes,Mc), ...
              'Color', color, 'LineStyle', lineStyle, 'LineWidth', 1.5, ...
              'DisplayName', legendEntry )

        if isLogScale
            set(gca,'YScale','log');
        end    

        h_leg = legend('show');
        set(h_leg,'Location','NorthWest','Interpreter','none')
        title(['Direct Aftershocks within ', num2str(nDays), ' days'])
        ylabel('Estimated aftershocks productivity')
        xlabel('Magnitude of triggering event (Mw)')
        
    end

end