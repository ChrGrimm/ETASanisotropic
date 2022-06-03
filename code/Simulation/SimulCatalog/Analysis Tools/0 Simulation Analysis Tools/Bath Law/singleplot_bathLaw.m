function singleplot_bathLaw( magnitudes, bathLaw_seq1, bathLaw_cat, bathLaw_historic, strTitle )

    % sequence results
    figure
    plot( magnitudes, bathLaw_seq1(:,2), 'b--', 'LineWidth', 1.5, 'DisplayName', 'Sequence' )
    hold on
%         plot( magnitudes, bathLaw_seq2(:,2), 'r--', 'LineWidth', 1.5, 'DisplayName', 'Full sequence' )
    fill( [magnitudes; flipud(magnitudes)], [bathLaw_seq1(:,1); flipud(bathLaw_seq1(:,3))], ...
          'blue', 'FaceAlpha', 0.2, 'EdgeColor', 'k', 'DisplayName', '5/95% (sequence)')

    % catalog results
    plot( magnitudes, bathLaw_cat(:,2), 'b', 'LineStyle', '-', 'LineWidth', 1.5, 'DisplayName', 'Catalog' )

    % Original catalog data
    scatter(bathLaw_historic(:,1), bathLaw_historic(:,2), 5*bathLaw_historic(:,3), 'k', 'filled', 'DisplayName', 'orig. catalog')

    % Ideal Bath Law line
    plot( magnitudes([1,end]), [1.2; 1.2], 'r:', 'LineWidth', 1.5, 'DisplayName', 'Bath`s Law')

    % Legend and title
    title(strTitle)
    legend show
    legend('Location', 'Northwest')
%     title('Bath Law')
    ylabel('Magnitude difference (Mw)')
    xlabel('Magnitude of triggering event (Mw)')
    ylim([0,4])
    
end