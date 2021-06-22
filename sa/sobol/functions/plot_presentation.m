function [] = plot_presentation(results, neurons, legend_cell, c_cell, output_array)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

set(0,'defaultAxesFontSize',40);

t_array = ["indices","ranking"];

%% Comparison figures
FIG = figure();

% for i = 1:2
%     subplot(2, 1, i);
    b = bar(c_cell{1});
    xticks(1:length(legend_cell)); 
    xticklabels(legend_cell); xlabel('Parameters', 'FontSize', 30); xtickangle(0)

    % Pretty colours
    colours = [colour_pairs('spring')];
    b(1).FaceColor = colours{1};
    b(2).FaceColor = colours{2};
    b(3).FaceColor = colours{3};

%     if i==1
        legend(output_array,'interpreter','latex', 'Location', 'best', 'FontSize', 35)
%     else
%         for j = 1:3
%             xtips = b(j).XEndPoints;
%             ytips = b(j).YEndPoints;
%             labels = string(b(j).YData);
%             text(xtips,ytips,labels,'HorizontalAlignment','center',...
%                 'VerticalAlignment','bottom', 'FontSize', 20)
%         end
%         yline(5,'b--')
%     end
% end

title(strcat("\textbf{Sobol FO} according to different model indicators"), 'FontSize', 45);
pause;
saveas(FIG,'figures/presentation/sobol_fo_comparison','epsc')

%% 

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

FIG = figure('Renderer', 'painters');

for i = 1:N
    state = results(:,idx_1(i), 1);
    p(i) = semilogx(neurons, state, '-s', ...
        'Color', colours_rnked(idx_1(i), :), ...
        'MarkerFaceColor', colours_rnked(idx_1(i), :), ...
        'LineWidth', 4, 'MarkerSize', 12); hold on;
    if i==1
        fill(CI_x, CI_y, 'black', 'linestyle', ':', 'facealpha', 0.05); 
        p(i) = semilogx(neurons, state, '-s', ...
            'Color', colours_rnked(idx_1(i), :), ...
            'MarkerFaceColor', colours_rnked(idx_1(i), :), ...
            'LineWidth', 4, 'MarkerSize', 12); hold on;    
    end
end
hold off;
ylim([min(0, lb),ub]); yticks(0:0.2:0.4);
xlim([neurons(1)-0.2, neurons(end)+50]); xticks([10 100 1000]); xlabel('Size of the network', 'FontSize', 30);
legend(p(1:2), legend_cell{idx_1(1:2)},'interpreter','latex', 'NumColumns', 1, 'FontSize', 50, 'Location', 'best');
title('\textbf{Sobol FO} convergence for $\gamma$','interpreter','latex', 'FontSize', 50)

pause;
exportgraphics(FIG,strcat('figures/presentation/sobol_fo_conv.eps'),'ContentType','vector','BackgroundColor','none')

set(0,'defaultAxesFontSize',35);
end

