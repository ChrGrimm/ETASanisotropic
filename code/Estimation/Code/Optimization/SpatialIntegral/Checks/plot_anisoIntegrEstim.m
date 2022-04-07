function plot_anisoIntegrEstim( iSample_integrAniso, ...
                                iSample_integrIso, ...
                                iSample_gradAniso_dD, ...
                                iSample_gradIso_dD, ...
                                iSample_gradAniso_dQ, ...
                                iSample_gradIso_dQ, ...
                                distGrid_integr, ...
                                distGrid_grad, ...
                                iMag, ...
                                iRupL, ...
                                iR )

%     D           = spatParamETAS(2);
%     gamma       = spatParamETAS(1);
%     q           = spatParamETAS(3);
%     magnitudes  = [4.85, 6.5, 8.0, 9.1]-4.85;
%     limitMesh   = [0.1, 0.4, 2, 6];
%     iRupLength   = get_iRupLength(magnitudes+4.85, typeRupL)/111;

    close all
    figure('units','normalized','outerposition',[0 0 1 1])
    xlimits = [0, iR]; % [0, 2.5*rupL]

    subplot(1,3,1)
    plot( distGrid_integr, iSample_integrAniso )
    hold on
    plot( distGrid_integr, iSample_integrIso )
    plot( distGrid_integr, iSample_integrAniso + iSample_integrIso )
    plot( [ 0.5*iRupL, 0.5*iRupL ], [0 1], 'r--' )
    plot( [ 1.0*iRupL, 1.0*iRupL ], [0 1], 'r--' )
    plot( [ 1.5*iRupL, 1.5*iRupL ], [0 1], 'r--' )
    plot( [ 2.0*iRupL, 2.0*iRupL ], [0 1], 'r--' )
    plot( [ 2.5*iRupL, 2.5*iRupL ], [0 1], 'r--' )
    title(['Spatial integral for Mw = ', num2str(iMag)])
%     xlim([0, 2.5*iRupL])
    xlim(xlimits)
    legend('aniso', 'iso', 'iso+aniso', 'X*RupL', 'Location', 'Southeast')

    %% Derivative by D
    subplot(1,3,2)
    plot( distGrid_grad, iSample_gradAniso_dD )
    hold on
    plot( distGrid_grad, iSample_gradIso_dD )
    plot( distGrid_grad, iSample_gradAniso_dD + iSample_gradIso_dD )
    plot( [ 0.5*iRupL, 0.5*iRupL ], [0 1], 'r--' )
    plot( [ 1.0*iRupL, 1.0*iRupL ], [0 1], 'r--' )
    plot( [ 1.5*iRupL, 1.5*iRupL ], [0 1], 'r--' )
    plot( [ 2.0*iRupL, 2.0*iRupL ], [0 1], 'r--' )
    plot( [ 2.5*iRupL, 2.5*iRupL ], [0 1], 'r--' )
    title(['Derivative of F by D for Mw = ', num2str(iMag)])
    xlim(xlimits)
    legend('aniso', 'iso', 'iso+aniso', 'X*RupL', 'Location', 'Southeast')

    %% Derivative by q
    subplot(1,3,3)
    plot( distGrid_grad, iSample_gradAniso_dQ )
    hold on
    plot( distGrid_grad, iSample_gradIso_dQ )
    plot( distGrid_grad, iSample_gradAniso_dQ + iSample_gradIso_dQ )
    plot( [ 0.5*iRupL, 0.5*iRupL ], [0 1], 'r--' )
    plot( [ 1.0*iRupL, 1.0*iRupL ], [0 1], 'r--' )
    plot( [ 1.5*iRupL, 1.5*iRupL ], [0 1], 'r--' )
    plot( [ 2.0*iRupL, 2.0*iRupL ], [0 1], 'r--' )
    plot( [ 2.5*iRupL, 2.5*iRupL ], [0 1], 'r--' )
    title(['Derivative of F by q for Mw = ', num2str(iMag)])
    xlim(xlimits)
    legend('aniso', 'iso', 'iso+aniso', 'X*RupL', 'Location', 'Southeast')

        %% Spatial density f
    %     figure(4)
    %     subplot(2,2,i)
    %     [X,Y] = meshgrid(-limitMesh(i):limitMesh(i)/100:limitMesh(i));
    %     r = get_dist2line([-iRupLength(i)/2, 0], [ iRupLength(i)/2, 0], [X(:),Y(:)]', 0)';
    %     f = (q-1)/(D*exp(gamma*magnitudes(i))) * f_inn(r, r.^2, iRupLength(i), magnitudes(i)).^(-q);
    %     f = reshape(f, size(X));
    %     h = surf(X, Y, f);
    % %         set(h,'LineStyle','none')
    %     title(['PDF spatial distribution for Mw = ', num2str(magnitudes(i)+4.85)])

end