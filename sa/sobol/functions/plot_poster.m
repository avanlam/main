function [] = plot_poster(results, neurons, legend_cell, output)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

N = length(legend_cell); colours_rnked = zeros(N, 3);

[~, idx_1] = sort(results(end, :, 1));

% Graph for INDICES : highest is most important
CI_x =[neurons, fliplr(neurons)]; 
i = 3; CI_y=[results(:,idx_1(i), 2)', flipud(results(:,idx_1(i), 3))'];

i1 = 12; i2 = 16;
colours_rnked(idx_1(1:i1), :) = colour_groups('blue', length(1:i1));
colours_rnked(idx_1(i1+1:i2), :) = colour_groups('green', length(i1+1:i2));
colours_rnked(idx_1(i2+1:end), :) = colour_groups('orange', length(i2+1:N));

idx_1 = flip(idx_1);

lb = min(results(:,:, 1), [], 'all');
ub = 0.4;

%% Visualisiation
set(0,'defaultAxesFontSize',45);

FIG = figure('Renderer', 'painters');
FIG.WindowState = 'fullscreen';

for i = 1:N
    state = results(:,idx_1(i), 1);
    p(i) = semilogx(neurons, state, '-s', ...
        'Color', colours_rnked(idx_1(i), :), ...
        'MarkerFaceColor', colours_rnked(idx_1(i), :), ...
        'LineWidth', 4, 'MarkerSize', 25); hold on;
    if i==1
        fill(CI_x, CI_y, 'black', 'linestyle', ':', 'facealpha', 0.05); 
        p(i) = semilogx(neurons, state, '-s', ...
            'Color', colours_rnked(idx_1(i), :), ...
            'MarkerFaceColor', colours_rnked(idx_1(i), :), ...
            'LineWidth', 4, 'MarkerSize', 25); hold on;    
    end
end
hold off;
ylim([min(0, lb),ub]); ylabel('Sobol FO', 'FontSize', 60); yticks(0:0.1:0.5);
xlim([neurons(1)-0.2, neurons(end)+50]); xticks([10 100 1000]); xlabel('Size of the network', 'FontSize', 60);
legend(p(1:2), legend_cell{idx_1(1:2)},'interpreter','latex', 'NumColumns', 1, 'FontSize', 60, 'Location', 'eastoutside');

exportgraphics(FIG,strcat('figures/', output, '/poster_fo.eps'),'ContentType','vector','BackgroundColor','none')
exportgraphics(FIG,strcat('figures/', output, '/poster_fo.png'))
% saveas(FIG,strcat('figures/', output, '/poster_fo'),'epsc')

set(0,'defaultAxesFontSize',35);
end

