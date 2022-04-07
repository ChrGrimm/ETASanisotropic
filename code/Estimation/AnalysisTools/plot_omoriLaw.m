function plot_omoriLaw( pathResultsFile, nDays, color, lineStyle, legendEntry )
% This function plots the estimated Omori Law, given an ETAS model results .mat file
%
% INPUTS:
% - pathResultsFile:        string      full path to .mat file containing ETAS estimation results
% - nDays:                  scalar      duration (in days) of aftershock triggering (e.g. 7, 100, 365)
% - color:                              line color code (Matlab notation)
% - lineStyle:                          line style code (Matlab notation)
% - legendEntry:            string      name of model to appear in legend

    %% Load or set ETAS parameters
    if ~isempty(pathResultsFile)
        % Load parameter estimates from submitted path to results file
        load(pathResultsFile, 'paramETAS')
        [~, ~, ~, c, p] = extract_paramETASvalues( paramETAS, 'default' );
    
    else
        % Manually set parameter estimates
        c       = 0.01;
        p       = 0.9;
    end
    
    %% Compute aftershock productivity function
    
    %% Plotting
    hold on
    t = 0:0.01:nDays;
    plot( t, (t+c).^-p, 'Color', color, 'LineStyle', lineStyle, 'LineWidth', 1.5, 'DisplayName', legendEntry )
    
    h_leg = legend('show');
    set(h_leg,'Location','NorthEast','Interpreter','none')
    title(['Omori Law (t+c)^{(-p)}'])
    ylabel('Aftershock Rate')
    xlabel('Time since mainshock (in days)')

end