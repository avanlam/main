function [FIG] = plot_error(x_all, y_all, xx, yy, lb, ub)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

colours = colour_pairs('spring');
    
%% Error

err_x_rom= 1/sqrt(numel(xx{1})).*reshape(vecnorm(vecnorm(xx{1}-xx{2},2, 1), 2, 2), size(x_all));
err_x_prom= 1/sqrt(numel(xx{1})).*reshape(vecnorm(vecnorm(xx{1}-xx{3},2, 1), 2, 2), size(x_all));
err_y_rom= 1/sqrt(numel(yy{1})).*reshape(vecnorm(vecnorm(yy{1}-yy{2},2, 1), 2, 2), size(y_all));
err_y_prom= 1/sqrt(numel(yy{1})).*reshape(vecnorm(vecnorm(yy{1}-yy{3},2, 1), 2, 2), size(y_all));
err_x_fom= zeros(size(err_x_prom)); err_y_fom= zeros(size(err_y_prom));
    
%% Graph
 
FIG = figure();

t = tiledlayout(1,2);
t.TileSpacing = 'tight';
t.TileIndexing = 'columnmajor';

    nexttile;
    ax1 = plot(x_all,err_x_fom, 'Color', colours{1}, 'Marker', 'x'); hold on; ax2 = plot(x_all,err_x_rom, 'Color', colours{2}, 'Marker', 'x'); ax3 = plot(x_all,err_x_prom, 'Color', colours{3}, 'Marker', 'x'); 
    xlabel('$\nu_6$'); ylabel('$|| err ||_2$'); xlim([lb(1) ub(1)]); ylim([0 0.02]);

    nexttile;
    plot(y_all,err_y_fom, 'Color', colours{1}, 'Marker', 'x'); hold on; plot(y_all,err_y_rom, 'Color', colours{2}, 'Marker', 'x'); plot(y_all,err_y_prom, 'Color', colours{3}, 'Marker', 'x'); 
    xlabel('$\Omega$'); set(gca,'ytick',[]); xlim([lb(2) ub(2)]); ylim([0 0.02]);

    lg = legend([ax1, ax2, ax3], ["FOM", "$ROM$","$pROM$"], 'interpreter', 'latex');
    lg.Layout.Tile = 'north';
    title(t, 'Robustness to variations in parameters', 'interpreter', 'latex', 'FontSize', 40)
    FIG.WindowState = 'fullscreen';
end

