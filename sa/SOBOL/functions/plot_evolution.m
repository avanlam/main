function [FIG] = plot_evolution(name, results, neurons, legendCell, tit, method, output)

N = length(legendCell); colours_rnked = zeros(N, 3);

[~, idx] = sort(results(end, :, 1));

FIG = figure('Renderer', 'painters');
FIG.WindowState = 'fullscreen';

if size(results, 3) > 1         % Graph for INDICES : highest is most important
    i = 3;
    CI_x =[neurons, fliplr(neurons)];
    CI_y=[results(:,idx(i), 2)', flipud(results(:,idx(i), 3))'];
    fill(CI_x, CI_y, 'black', 'linestyle', ':', 'facealpha', 0.05); hold on;
else                            % Graph for RANKING : lowest is most important
    idx = flip(idx);
end

colours_rnked(idx(1:floor(2*N/3)), :) = colour_groups('blue', length(1:floor(2*N/3)));
colours_rnked(idx(ceil(2*N/3):floor(5*N/6)), :) = colour_groups('green', length(ceil(2*N/3):floor(5*N/6)));
colours_rnked(idx(ceil(5*N/6):end), :) = colour_groups('orange', length(ceil(5*N/6):N));

if size(results, 3) > 1
    idx = flip(idx);
end


for i = 1:N
    state = results(:,idx(i), 1);
    p(i) = semilogx(neurons, state, '-s', 'Color', colours_rnked(idx(i), :)); hold on;
end
hold off;
title(['\textbf{Convergence of the Sobol ', tit, '}'], 'FontSize', 45);
subtitle(['according to the indicator \verb|', output, '|'], 'FontSize', 45, 'interpreter', 'latex');
xlabel("Size of the neuronal network", 'FontSize', 40); ylabel(['Sobol ', tit], 'FontSize', 40);
legend(p, legendCell{idx},'interpreter','latex', 'NumColumns', 1, 'Location', 'best', 'FontSize', 30)

saveas(FIG,strcat('figures_', method, '/', output, '/', name),'epsc')

% FIG = figure();
% FIG.WindowState = 'fullscreen';
% b = bar3(results(:,idx)');
% title(['\textbf{Convergence of the Sobol ', tit, '}'], 'FontSize', 45);
% subtitle(['according to the indicator \verb|', output, '|'], 'FontSize', 45, 'interpreter', 'latex');
% xticks(1:length(neurons)); xticklabels(neurons); xtickangle(90)
% xlabel("Size", 'FontSize', 40, 'rotation', 30); 
% ytickangle(90); ylim([-1, length(legendCell)+1]);
% yticks(1:length(legendCell)); yticklabels(legendCell(idx));
% zlabel(['Sobol ',tit], 'FontSize', 40); zticks(0:2:16);
% colormap spring %summer
% 
% for k = 1:length(b)
%     zdata = b(k).ZData;
%     b(k).CData = zdata;
%     b(k).FaceColor = 'interp';
% end
% 
% saveas(FIG,strcat(method, '/', output, '/', name, '_3D'),'epsc')

end

