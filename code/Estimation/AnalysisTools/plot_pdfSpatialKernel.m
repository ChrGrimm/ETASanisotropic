function plot_pdfSpatialKernel( pathResultsFile, magnitudes )
% This function plots the 3D probability density function of a spatial kernel, given an ETAS model results .mat file.
% Spatial X- and Y-axis are in normalized to units of rupture lengths.
%
% INPUTS:
% - pathResultsFile:        string              full path to .mat file containing ETAS estimation results
% - magnitudes:             numerical vector    vector of magnitudes, for which pdfs are plotted in respective subplots
    
    %% Load ETAS results
    % Load parameter estimates from submitted path to results file
    load(pathResultsFile, 'paramETAS', 'EstimSettings', 'ModelFuncs')
    [~,~,~,~,~,~,~,q]   = extract_paramETASvalues( paramETAS, 'default' );
    anisoFromMw         = EstimSettings.SpaceSettings.anisoFromMw;
    restrFactor         = EstimSettings.SpaceSettings.restrFactor;
    
    %% Loop over magnitudes
    nMagnitudes = length(magnitudes);
    for iMag = 1:nMagnitudes
        
        subplot(1,nMagnitudes,iMag)

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
    
        % Define 2D spatial mesh grid
        meshExtent  = min(3,restrFactor) * rupL;
        [X,Y]       = meshgrid(-meshExtent:meshExtent/50:meshExtent);
        
        % Compute distances to rupture
        if strcmp(typeKernel, 'aniso')
            r = get_dist2line([-rupL/2, 0], [ rupL/2, 0], [X(:),Y(:)]', 'degree', 0)';

        else
            r = getDistance(X(:), Y(:), 0, 0, 'degree');

        end
    
        % Compute pdf on grid X, Y
        f_pdf           = (r <= spatRestr) .* ModelFuncs.fSpatK_factor(spatRestr, Mw, rupExtent) ...
                            .* ModelFuncs.fSpatK_inn(r, r.^2, Mw, rupExtent).^(-q);
        f_pdf           = reshape(f_pdf, size(X));
        f_pdf(f_pdf==0) = NaN;

        %% Plotting
        surf(X/rupL, Y/rupL, f_pdf);

        %% Figure design
        zlim([0, max(max(f_pdf))])
        axisLimit = min(3,restrFactor);
        xlim([-axisLimit, axisLimit])
        ylim([-axisLimit, axisLimit])

        xlabel('Normalized space (x)')
        ylabel('Normalized space (y)')
        zlabel('pdf')

    %     make_titleLeftCorner( '(a)' )
        title(['Mw = ', num2str(Mw)])

    %     a = get(gca,'ZTickLabel');  
    %     set(gca,'ZTickLabel',a,'fontsize',nFontSize)
    
    end
end