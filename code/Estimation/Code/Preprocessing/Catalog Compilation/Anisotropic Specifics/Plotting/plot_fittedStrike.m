function plot_fittedStrike( isInTime, ...
                              idxAniso, ...
                              strikes, ...
                              epiPos, ...
                              x, y, t, mag, wi, rupL, ...
                              lon, lat, ...
                              spaceUnit, ...
                              spatWindow, ...
                              minSpaceDist, ...
                              g_inn, f_inn, p, q )
                                        
    iAnisoEv = 0;
    for iEvent = idxAniso'

        iAnisoEv                = iAnisoEv + 1; 
        forwardRateContrib_sum  = cell(1,length(strikes{iAnisoEv}));
        forwardRateContrib_logl = cell(1,length(strikes{iAnisoEv}));
        f_inn_ij                = cell(1,length(strikes{iAnisoEv}));
        g_inn_ij                = cell(1,length(strikes{iAnisoEv}));
        lambda_trig             = cell(1,length(strikes{iAnisoEv}));
        jEvTimeSpace            = cell(1,length(strikes{iAnisoEv}));
        dist                    = cell(1,length(strikes{iAnisoEv}));
        
        for iStrike = 1:length(strikes{iAnisoEv})
            
            % Compute start and end points of rupture segments
            [ rupStart(iEvent,:), ...
                rupEnd(iEvent,:) ] = compute_rupCoords( x( iEvent ), ...
                                                      y( iEvent ), ...
                                                      rupL( iEvent ), ...
                                                      strikes{iAnisoEv}(iStrike), ...
                                                      epiPos{iAnisoEv}(iStrike) );

            [ forwardRateContrib_sum{iStrike}, ...
               forwardRateContrib_logl{iStrike}, ...
               ~, ...
               f_inn_ij{iStrike}, ...
               g_inn_ij{iStrike}, ....
               lambda_trig{iStrike}, ...
               jEvTimeSpace{iStrike}, ...
               dist{iStrike} ] = estimate_forwardRateContrib( isInTime, ...
                                                              iAnisoEv, ...
                                                              x, y, t, mag, wi, rupL, ...
                                                              iEvent, ...
                                                              spaceUnit, ...
                                                              spatWindow, ...
                                                              rupStart, rupEnd, ...
                                                              minSpaceDist, ...
                                                              g_inn, f_inn, p, q );
                                                          
        end
        
        forwardRateContrib_sum  = 1/length(strikes{iAnisoEv}) * sum(cell2mat(forwardRateContrib_logl),2);
        forwardRateContrib_logl = 1/length(strikes{iAnisoEv}) * sum(cell2mat(forwardRateContrib_logl),2);
        f_inn_ij                = 1/length(strikes{iAnisoEv}) * sum(cell2mat(f_inn_ij),2);
        g_inn_ij                = 1/length(strikes{iAnisoEv}) * sum(cell2mat(g_inn_ij),2);
        lambda_trig             = 1/length(strikes{iAnisoEv}) * sum(cell2mat(lambda_trig),2);
        jEvTimeSpace            = jEvTimeSpace{1};
        distKm                  = 111 * min(cell2mat(dist),[],2);

        figure
        lambda_plot = max(5, 50*lambda_trig); %(log10(lambda_trig)+10).^2; %.^2 10*f_inn_ij.^(-q)
%         scatter( lon(jEvTimeSpace), lat(jEvTimeSpace), lambda_plot, log10(g_inn_ij), 'filled' )  % (2./distKm).^2
        scatter( lon(jEvTimeSpace), lat(jEvTimeSpace), 10, log10(g_inn_ij), 'filled' )
%         scatter( x(jEvTimeSpace), y(jEvTimeSpace), lambda_plot, g_inn_ij.^(-q), 'filled' )
        colorbar
        
        for iStrike = 1:length(strikes{iAnisoEv})
            thisStrike = strikes{iAnisoEv}(iStrike);
            iEpiPos = epiPos{iAnisoEv}(iStrike);
            [ rupS_plot, rupE_plot ] = compute_rupCoords( lon(iEvent), lat(iEvent), rupL(iEvent), thisStrike, iEpiPos );
            xy_low = [lon(iEvent), lat(iEvent)] - minSpaceDist .* [ cosd(180-thisStrike), sind(180-thisStrike) ];
            [ rupS_low, rupE_low ] = compute_rupCoords( xy_low(1), xy_low(2), rupL(iEvent), thisStrike, iEpiPos );
            xy_upp = [lon(iEvent), lat(iEvent)] + minSpaceDist .* [ cosd(180-thisStrike), sind(180-thisStrike) ];
            [ rupS_upp, rupE_upp ] = compute_rupCoords( xy_upp(1), xy_upp(2), rupL(iEvent), thisStrike, iEpiPos );
            hold on 
            plot( lon(iEvent), lat(iEvent), 'r*', 'LineWidth', 10)
%             plot( xy_low(1), xy_low(2), 'r*' )
%             plot( xy_upp(1), xy_upp(2), 'r*' )
            plot([rupS_plot(1),rupE_plot(1)], [rupS_plot(2),rupE_plot(2)], 'r-', 'LineWidth', 2)
            plot([rupS_low(1),rupE_low(1)], [rupS_low(2),rupE_low(2)], 'r--')
            plot([rupS_upp(1),rupE_upp(1)], [rupS_upp(2),rupE_upp(2)], 'r--')
        end
        axis equal
        xlabel('Longitude')
        ylabel('Latitude')
        make_titleLeftCorner( '(c)' )
%         title(['Strike = ', num2str(thisStrike), ', EpiPos = ', num2str(iEpiPos), ...
%                ', Sum = ', num2str(forwardRateContrib_sum), ', Log-L = ', num2str(forwardRateContrib_logl)])
%         [ sort_lambda, idxSort ] = sort(lambda_trig, 'descend');
%         sort_ginn                = g_inn_ij( idxSort );
%         sort_finn                = sort( f_inn_ij.^(-q) );
%         test                     = [lambda_plot, g_inn_ij.^(-q), f_inn_ij.^(-p)];
        
    end

end

%% Old code pieces
    %     plot(Catalog.lon(idxSortedRates), Catalog.lat(idxSortedRates), 'r*', 'LineWidth', 5)
    %     axis([floor(min(polygonComple(:,1))), ceil(max(polygonComple(:,1))), ...
    %           floor(min(polygonComple(:,2))), ceil(max(polygonComple(:,2)))])
    %     plot(polygonTarget(:,1), polygonTarget(:,2), 'r-')
    %     plot(polygonComple(:,1), polygonComple(:,2), 'r--')
