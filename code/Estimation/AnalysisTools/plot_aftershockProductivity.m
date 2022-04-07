function plot_aftershockProductivity( pathResultsFile, nDays, color, lineStyle, legendEntry, isLogScale )
% This function plots the aftershock productivity, given an ETAS model results .mat file
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
        load(pathResultsFile, 'paramETAS', 'EstimSettings')
        [~, A, alpha, c, p] = extract_paramETASvalues( paramETAS, 'default' );
        Mc                  = EstimSettings.TargetSettings.Mc;
        Mmax                = EstimSettings.TargetSettings.Mmax;
    
    else
        % Manually set parameter estimates
        A       = 0.01;
        alpha   = 2;
        c       = 0.01;
        p       = 0.9;
        Mc      = 4;
        Mmax    = 9;
    end
    
    %% Compute aftershock productivity function
    % As the Omori law is not normalized (does not integrate to 1), scale aftershock productivity
    OmoriIntegral       = 1/(1-p) * ( (nDays+c).^(1-p) - c.^(1-p) );
    funcProductivity    = @(mag, Mc) OmoriIntegral * A * exp(alpha*(mag-Mc));
    % Magnitude vector
    magnitudes = Mc:0.1:Mmax;
    
    %% Plotting
    hold on
    plot( magnitudes, funcProductivity(magnitudes,Mc), ...
          'Color', color, 'LineStyle', lineStyle, 'LineWidth', 1.5, ...
          'DisplayName', legendEntry )
      
    if isLogScale
        set(gca,'YScale','log');
    end
    
    h_leg = legend('show');
    set(h_leg,'Location','NorthWest','Interpreter','none')
    title(['Direct Aftershocks within ', num2str(nDays), ' days'])
    ylabel('Estimated aftershock productivity')
    xlabel('Magnitude of triggering event (Mw)')

end