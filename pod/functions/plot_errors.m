function [FIG] = plot_errors(name, x_all, y_all, z_all, xx, yy, zz, xx_omega, yy_omega, zz_omega, lb, ub)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

colours = colour_pairs('springXL');
    
%% Error pROM

err_x_rom= 1/sqrt(numel(xx{1})).*reshape(vecnorm(vecnorm(xx{1}-xx{2},2, 1), 2, 2), size(x_all));
err_x_prom= 1/sqrt(numel(xx{1})).*reshape(vecnorm(vecnorm(xx{1}-xx{3},2, 1), 2, 2), size(x_all));

err_y_rom= 1/sqrt(numel(yy{1})).*reshape(vecnorm(vecnorm(yy{1}-yy{2},2, 1), 2, 2), size(y_all));
err_y_prom= 1/sqrt(numel(yy{1})).*reshape(vecnorm(vecnorm(yy{1}-yy{3},2, 1), 2, 2), size(y_all));

err_z_rom= 1/sqrt(numel(zz{1})).*reshape(vecnorm(vecnorm(zz{1}-zz{2},2, 1), 2, 2), size(z_all));
err_z_prom= 1/sqrt(numel(zz{1})).*reshape(vecnorm(vecnorm(zz{1}-zz{3},2, 1), 2, 2), size(z_all));

err_x_fom= zeros(size(err_x_prom)); err_y_fom= zeros(size(err_y_prom)); err_z_fom= zeros(size(err_z_prom)); 
    
%% Error omega pROM

err_x_prom_omega = 1/sqrt(numel(xx_omega{1})).*reshape(vecnorm(vecnorm(xx_omega{1}-xx_omega{3},2, 1), 2, 2), size(x_all));
err_y_prom_omega= 1/sqrt(numel(yy_omega{1})).*reshape(vecnorm(vecnorm(yy_omega{1}-yy_omega{3},2, 1), 2, 2), size(y_all));
err_z_prom_omega= 1/sqrt(numel(zz_omega{1})).*reshape(vecnorm(vecnorm(zz_omega{1}-zz_omega{3},2, 1), 2, 2), size(z_all));
    
%% Graph
 
FIG = figure();

t = tiledlayout(1,3);
t.TileSpacing = 'tight';

    % pROM
    nexttile;
    ax1 = plot(x_all,err_x_fom, 'Color', colours{1}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); hold on; 
    ax2 = plot(x_all,err_x_rom, 'Color', colours{2}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); 
    ax3 = plot(x_all,err_x_prom, 'Color', colours{3}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); 
    ax4 = plot(x_all,err_x_prom_omega, 'Color', colours{4}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); 
    xlabel('$\Omega$'); ylabel('$|| err ||_2$'); xlim([lb(1) ub(1)]); % ylim([0 0.02]);

    nexttile;
    plot(y_all,err_y_fom, 'Color', colours{1}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); hold on; 
    plot(y_all,err_y_rom, 'Color', colours{2}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); 
    plot(y_all,err_y_prom, 'Color', colours{3}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); 
    plot(y_all,err_y_prom_omega, 'Color', colours{4}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); 
    xlabel('$L_0$'); set(gca,'ytick',[]); xlim([lb(2) ub(2)]); % ylim([0 0.02]);

    nexttile;
    plot(z_all,err_z_fom, 'Color', colours{1}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); hold on; 
    plot(z_all,err_z_rom, 'Color', colours{2}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); 
    plot(z_all,err_z_prom, 'Color', colours{3}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); 
    plot(z_all,err_z_prom_omega, 'Color', colours{4}, 'LineWidth', 3, 'Marker', 'x', 'MarkerSize', 18); 
    xlabel('$\nu_6$'); set(gca,'ytick',[]); xlim([lb(3) ub(3)]); % ylim([0 0.02]);

    lg = legend([ax3, ax4, ax1, ax2], ["$pROM$", name, "FOM", "$ROM$"], 'interpreter', 'latex', 'NumColumns', 2);
    lg.Layout.Tile = 'north';
    title(t, '\textbf{Robustness to variations in parameters}', 'interpreter', 'latex', 'FontSize', 50)
    FIG.WindowState = 'fullscreen';
end

