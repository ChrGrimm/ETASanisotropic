function sharedplot_bathLaw( magnitudes, bathLaw_seq1, bathLaw_cat, bathLaw_historic, strLegend, strRegion, iRun )

    color       = define_colors4plotsInLoops( iRun );
    % sequence results
    plot( magnitudes, bathLaw_seq1(:,2), 'Color', color, 'LineStyle', '--', ...
        'LineWidth', 1.5, 'DisplayName', [strLegend, ' synth. sequences'] )
    hold on
%         plot( magnitudes, bathLaw_seq2(:,2), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Full sequence' )
%         if iRun == 1
%             fill( [magnitudes; flipud(magnitudes)], [bathLaw_seq1(:,1); flipud(bathLaw_seq1(:,3))], ...
%                 'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'DisplayName', [strLegend, ' sequence 5/95%'])
%         end

    % catalog results
    plot( magnitudes, bathLaw_cat(:,2), 'Color', color, 'LineStyle', '-', ...
        'LineWidth', 1.5, 'DisplayName', [strLegend, ' synth. catalogs'] )
    hold on

    if iRun == 2 % 2
%             fill( [magnitudes; flipud(magnitudes)], [bathLaw_seq1(:,1); flipud(bathLaw_seq1(:,3))], ...
%                 'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'DisplayName', [strLegend, ' sequence 5/95%'])
        fill( [magnitudes; flipud(magnitudes)], [bathLaw_cat(:,1); flipud(bathLaw_cat(:,3))], ...
            'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'DisplayName', [strLegend, ' synth. catalog 10/90%'])
    end

    if iRun == 2 % 2
        % Original catalog data
        if contains(strLegend, 'CAL')
            pointSize = 10;
        else
            pointSize = 5;
        end
        scatter(bathLaw_historic(:,1), bathLaw_historic(:,2), pointSize*bathLaw_historic(:,3), ...
                'k', 'filled', 'DisplayName', 'original catalog')
    end

    % Legend & Format
    legend show
    legend('Location', 'Northwest')

%     title('Bath Law')
    ylabel('Magnitude difference (Mw)')
    xlabel('Magnitude of triggering event (Mw)')
    if contains(strRegion, 'CAL')
        ylim([0,4.5])
        xlim([5.5, 7.5])
        text(6.9, 4, 'Southern California', 'FontSize', 10)
    else
        ylim([0,4.5])
        xlim([5.5, 9.0])
        text(8.5, 4, 'Japan', 'FontSize', 10)
    end
    
end