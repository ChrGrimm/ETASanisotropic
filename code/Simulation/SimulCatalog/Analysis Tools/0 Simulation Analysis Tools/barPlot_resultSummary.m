close all
figure
bar(table2array(nEvents))
ylim([0,1])
hold on
plot([0.5, 3.5], [0.025, 0.025], 'k--', 'LineWidth', 1.5)
plot([0.5, 3.5], [0.975, 0.975], 'k--', 'LineWidth', 1.5, 'DisplayName', 'none')
legend({'ETAS convent.', 'ETAS iso-r', 'ETAS aniso-r', 'ETASI iso-r', 'ETASI aniso-r', '2.5% / 97.5%'}, ...
    'Location', 'Northwest')
xticklabels({'Experiment 1', 'Experiment 2', 'Experiment 3'})
ylabel('P(#Events >= observation)')
make_titleLeftCorner('(a)')
set(findall(gcf,'-property','FontSize'),'FontSize',13)

figure
bar(table2array(largestAftershock))
ylim([0,1])
hold on
plot([0.5, 3.5], [0.025, 0.025], 'k--', 'LineWidth', 1.5)
plot([0.5, 3.5], [0.975, 0.975], 'k--', 'LineWidth', 1.5, 'DisplayName', 'none')
% legend({'ETAS conventional', 'ETAS iso-r', 'ETAS aniso-r', 'ETASI iso-r', 'ETASI aniso-r', '2.5%/97.5% Quantiles'}, ...
%     'Location', 'Northwest')
xticklabels({'Experiment 1', 'Experiment 2', 'Experiment 3'})
ylabel('P(largest aftershock >= observation)')
make_titleLeftCorner('(b)')
set(findall(gcf,'-property','FontSize'),'FontSize',13)

figure
bar(-1,0.5, 'DisplayName', 'none')
hold on
b = bar(table2array(InfoGains));
xlim([0.5,3.5])
legend(b, {'ETAS iso-r', 'ETAS aniso-r', 'ETASI iso-r', 'ETASI aniso-r'}, ...
    'Location', 'Northwest')
xticks([1,2,3])
xticklabels({'Experiment 1', 'Experiment 2', 'Experiment 3'})
ylabel('Information Gain Relative to ETAS conv.')
make_titleLeftCorner('(c)')
set(findall(gcf,'-property','FontSize'),'FontSize',13)

% figure
% b = bar(3 + [-0.2, -0.1, 0, 0.1, 0.2], table2array(InfoGains(3,:)));
% hold on
% bar(1, table2array(InfoGains(1,2:end)));
% bar(2, table2array(InfoGains(2,2:end)));
% b = bar(table2array(InfoGains));
% xlim([0.5,3.5])
% ylim([-1, 2])
% legend(b, {'ETAS conventional', 'ETAS iso-r', 'ETAS aniso-r', 'ETASI iso-r', 'ETASI aniso-r'}, ...
%     'Location', 'Northwest')
% xticks([1,2,3])
% xticklabels({'Experiment 1', 'Experiment 2', 'Experiment 3'})
% ylabel('Inform. Gain (relative to ETAS conv.)')
% make_titleLeftCorner('(c)')
% set(findall(gcf,'-property','FontSize'),'FontSize',13)
