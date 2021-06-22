function [] = plot_poster_control(L0_sample, lossfcn_L0, iter_fom, iter_rom, iter_prom, delta_loss_fom, delta_loss_rom, delta_loss_prom, n_latent_pod)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

poster = 0;
if poster 
    colours{1} = [255,0,0]./255;
    colours{2} = [255,153,51]./255;
    colours{3} = [102,178,255]./255;
    c = [0.8 0.8 0.8];
    set(0,'defaultAxesFontSize',45);
else
    colours = colour_pairs('spring');
    c = 'k';
    set(0,'defaultAxesFontSize',40);
end
    
FIG = figure;

[~, idx_rom] = min(delta_loss_rom(:,2)); [~, idx_prom] = min(delta_loss_prom(:,2)); 

plot(L0_sample,lossfcn_L0, 'Color', c, 'LineWidth', 4); hold on;
ax1 = plot(iter_fom(:,2),delta_loss_fom(:,2),'x','color',colours{1}, 'MarkerSize', 25);
ax2 = plot(iter_rom(:,2),delta_loss_rom(:,2),'.','color',colours{2}, 'MarkerSize', 50);
ax3 = plot(iter_prom(:,2),delta_loss_prom(:,2),'.','color',colours{3}, 'MarkerSize', 50);
xline(iter_fom(end,2), 'color',colours{1}, 'LineWidth', 6); 
xline(iter_rom(idx_rom,2), '--', 'color',colours{2}, 'LineWidth', 6); 
xline(iter_prom(idx_prom,2), ':', 'color',colours{3}, 'LineWidth', 8);
xlim([0.005 0.02]);
xlabel('$L_0$','interpreter','latex');

if poster 
    ylabel('$\mathcal{L}(L_0)$'); 
else
    title('$\mathcal{L}(L_0)$','interpreter','latex', 'FontSize', 45); 
end

legend([ax1, ax2, ax3], {'FOM','ROM','pROM'},'interpreter','latex', 'FontSize', 40);
% title('\textbf{Optimisation: loss evolution}', 'interpreter','latex', 'FontSize', 60)

pause;
if poster 
    name = 'figures/poster_control.eps';
else
    name = 'figures/presentation/control.eps';
end
exportgraphics(FIG,strcat(name),'ContentType','vector','BackgroundColor','none')

set(0,'defaultAxesFontSize',35);
end