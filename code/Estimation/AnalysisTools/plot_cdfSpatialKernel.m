function plot_cdfSpatialKernel( pathResultsFile, magnitudes, color, linestyle, legendEntry )
% This function plots the cumulative density function of a spatial kernel against the distance (normalized by rupture length)
% to an isotropic or anisotropic rupture, given an ETAS model results .mat file
%
% INPUTS:
% - pathResultsFile:        string              full path to .mat file containing ETAS estimation results
% - magnitudes:             numerical vector    vector of magnitudes, for which cdfs are plotted in respective subplots
% - color:                                      line color code (Matlab notation)
% - lineStyle:                                  line style code (Matlab notation)
% - legendEntry:            string              name of model to appear in legend

    %% Load ETAS results
    % Load parameter estimates from submitted path to results file
    load(pathResultsFile, 'paramETAS', 'EstimSettings', 'ModelFuncs')
    [~,~,~,~,~,~,~,q]   = extract_paramETASvalues( paramETAS, 'default' );
    anisoFromMw         = EstimSettings.SpaceSettings.anisoFromMw;
    
    %% Loop over magnitudes
    nMagnitudes = length(magnitudes);
    for iMag = 1:nMagnitudes
        
        subplot(1,nMagnitudes,iMag)
        hold on

        %% Evaluate cdf of spatial kernel on distance grid
        % Compute normalization factor due to spatial restriction
        Mw = magnitudes(iMag);
        if Mw >= anisoFromMw
            typeKernel  = 'aniso';
        else
            typeKernel  = 'iso';
        end
        spatRestr       = ModelFuncs.fRestrSpace_degrees( Mw, typeKernel );
        rupL            = ModelFuncs.fRupLength( Mw );
        rupExtent       = rupL * strcmp(typeKernel, 'aniso');
        normalizFactor  = (1-ModelFuncs.fSpatK_inn(spatRestr, spatRestr.^2, Mw, rupExtent).^(1-q));
    
        % Evaluate cdf on distance grid
        r               = 0:rupL/1000:spatRestr;
        spatKernel_cdf  = (1 - ModelFuncs.fSpatK_inn(r, r.^2, Mw, rupExtent).^(1-q)) / normalizFactor;
    
        %% Plotting
        semilogx(r/rupL, spatKernel_cdf, 'Color', color, 'LineStyle', linestyle, 'LineWidth', 1.5, 'DisplayName', legendEntry)

        % Figure design
        h_leg = legend('show');
        set(h_leg,'Location','Southeast','Interpreter','none')
%         make_titleLeftCorner( '(a)' )
        title('Spatial kernel cdf by distance from rupture')
        xlabel('Distance (normalized by rupture length)')
        ylabel('Spatial cdf') % 'CDF'

    %     a = get(gca,'XTickLabel');  
    %     set(gca,'XTickLabel',a,'fontsize',nFontsize)
    
    end
end