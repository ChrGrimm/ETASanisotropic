function plot_triggRelations_ETASvsZBZ( TriggerRelationsETAS, TriggerRelationsZBZ )
% Compare ETAS and ZBZ trigger results

    % Comparison ZBZ vs. ETAS
    figure
    scatter(log10(TriggerRelationsZBZ.tij_rescaled), log10(TriggerRelationsZBZ.rij_rescaled), 5, ...
            TriggerRelationsETAS.backgrProb, 'filled')
    hold on
    xlims = [min(log10(TriggerRelationsZBZ.tij_rescaled)), max(log10(TriggerRelationsZBZ.tij_rescaled))];
    plot(xlims, -(xlims+3), 'r--')
    colorbar

end