function FIG1 = plot_optimal_projection(full, U_pod, V_pod)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

num_osc = size(U_pod,1)./4;

rom.x = (U_pod * U_pod' * full.x')';
prom.x = (V_pod * V_pod' * full.x')';
colours = colour_pairs('spring');


    FIG1 = figure();
    t = tiledlayout(2,2);
    t.TileSpacing = 'tight';
    sgtitle('\textbf{Oscillators in each model}', 'FontSize', 40, 'interpreter', 'latex');
    
    nexttile;
    ax1 = plot(1:num_osc,full.x(end,1:4:end-3),'.','color',colours{1}, 'MarkerSize', 30); ylabel({'$(X_i)$'},'interpreter','latex'); hold on;
    ax2 = plot(1:num_osc,rom.x(end,1:4:end-3), 'x', 'color',colours{2}, 'MarkerSize', 10);
    ax3 = plot(1:num_osc,prom.x(end,1:4:end-3), 'x','color',colours{3}, 'MarkerSize', 10); set(gca,'xtick',[])

    nexttile;
    plot(1:num_osc,full.x(end,2:4:end-2),'.','color',colours{1}, 'MarkerSize', 30); ylabel({'$(Y_i)$'},'interpreter','latex'); hold on;
    plot(1:num_osc,rom.x(end,2:4:end-2), 'x','color',colours{2}, 'MarkerSize', 10); 
    plot(1:num_osc,prom.x(end,2:4:end-2), 'x','color',colours{3}, 'MarkerSize', 10); set(gca,'xtick',[])
    
    nexttile;
    plot(1:num_osc,full.x(end,3:4:end-1),'.','color',colours{1}, 'MarkerSize', 30); ylabel({'$(Z_i)$'},'interpreter','latex'); hold on;
    plot(1:num_osc,rom.x(end,3:4:end-1), 'x','color',colours{2}, 'MarkerSize', 10);  
    plot(1:num_osc,prom.x(end,3:4:end-1), 'x','color',colours{3}, 'MarkerSize', 10); xlabel({'$i$'},'interpreter','latex');
    
    nexttile;
    plot(1:num_osc,full.x(end,4:4:end),'.','color',colours{1}, 'MarkerSize', 30); ylabel({'$(V_i)$'},'interpreter','latex'); hold on; 
    plot(1:num_osc,rom.x(end,4:4:end), 'x','color',colours{2}, 'MarkerSize', 10);
    plot(1:num_osc,prom.x(end,4:4:end), 'x','color',colours{3}, 'MarkerSize', 10); xlabel({'$i$'},'interpreter','latex'); 

    lg = legend([ax1(1), ax2(1), ax3(1)], ["FOM", "$ROM$","$pROM$"], 'interpreter', 'latex');
    lg.Layout.Tile = 'north';
end

