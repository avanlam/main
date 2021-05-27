function [FIG] = plot_evolution_together(name, results_1, results_2, neurons, legend_cell, tit_cell, output)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = length(legend_cell); colours_rnked = zeros(N, 3);

[~, idx_1] = sort(results_1(end, :, 1));

% Graph for INDICES : highest is most important
CI_x =[neurons, fliplr(neurons)]; 
i = 3; CI_y=[results_1(:,idx_1(i), 2)', flipud(results_1(:,idx_1(i), 3))'];

i1 = 5; i2 = 8;
colours_rnked(idx_1(1:i1), :) = colour_groups('blue', length(1:i1));
colours_rnked(idx_1(i1+1:i2), :) = colour_groups('green', length(i1+1:i2));
colours_rnked(idx_1(i2+1:end), :) = colour_groups('orange', length(i2+1:N));

idx_1 = flip(idx_1);

% Graph for RANKING : lowest is most important
idx_2 = idx_1;

fo = tit_cell{1}(1:2);

if strcmp( fo, 'FO' )
    lb = 1;
else
    lb = max(results_1(:,:, 1), [], 'all');
end
ub = min(results_1(:,:, 1), [], 'all');

%% Visualisiation
set(0,'defaultAxesFontSize',30);

FIG = figure('Renderer', 'painters');
FIG.WindowState = 'fullscreen';

t = tiledlayout(2,1);
t.TileSpacing = 'tight';

title(t, ['\textbf{Sobol ', fo,' } according to the indicator \verb|', output, '|'], 'FontSize', 50, 'interpreter', 'latex');

i1 = nexttile;
for i = 1:N
    state = results_1(:,idx_1(i), 1);
    p(i) = semilogx(neurons, state, '-s', 'Color', colours_rnked(idx_1(i), :)); hold on;
    if i==1
        fill(CI_x, CI_y, 'black', 'linestyle', ':', 'facealpha', 0.05); 
        p(i) = semilogx(neurons, state, '-s', 'Color', colours_rnked(idx_1(i), :));
    end
end
hold off;
ylim([min(0, ub),lb]); ylabel(tit_cell{1}, 'FontSize', 40);
set(gca,'xtick',[])

lg = legend(p, legend_cell{idx_1},'interpreter','latex', 'NumColumns', 1, 'FontSize', 35);
lg.Layout.Tile = 'east';

i2 = nexttile;
for i = 1:N
    state = results_2(:,idx_2(i), 1);
    p(i) = semilogx(neurons, state, '-s', 'Color', colours_rnked(idx_2(i), :)); hold on;
end
hold off;
ylim([0,17]);yticks(1:3:16);
xlabel("Size of the neuronal network", 'FontSize', 40); ylabel( tit_cell{2}, 'FontSize', 40);

linkaxes([i1 i2],'x')


saveas(FIG,strcat('figures/', output, '/', name),'epsc')
end

