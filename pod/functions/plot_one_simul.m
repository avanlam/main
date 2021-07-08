function [FIG1, FIG2] = plot_one_simul(num_n, n_asympt, full, rom, prom)
%PLOT_ONE_SIMUL draws the simulation solution and the synchronisation
%outputs of a single time simulation of the modified Goodwin model
    
if numel(n_asympt) == 1
    n_asympt = [n_asympt n_asympt n_asympt];
end

colours = colour_pairs('spring');

    Noscil_plot=floor(num_n/5);

    FIG1 = figure();
    t = tiledlayout(4,3);
    t.TileSpacing = 'tight';
    sgtitle('\textbf{Oscillators in each model}', 'FontSize', 50, 'interpreter', 'latex');
    
    nexttile(1);
    ax1 = plot(full.time,full.x(:,1:4*Noscil_plot:end-3),'color',colours{1}); ylabel({'$(X_i)_{FOM}$'},'interpreter','latex'); xlim([0 500]); set(gca,'xtick',[]) 
    nexttile(4);
    plot(full.time,full.x(:,2:4*Noscil_plot:end-2),'color',colours{1}); ylabel({'$(Y_i)_{FOM}$'},'interpreter','latex'); xlim([0 500]); set(gca,'xtick',[])
    nexttile(7);
    plot(full.time,full.x(:,3:4*Noscil_plot:end-1),'color',colours{1}); ylabel({'$(Z_i)_{FOM}$'},'interpreter','latex'); xlim([0 500]); set(gca,'xtick',[])
    nexttile(10);
    plot(full.time,full.x(:,4:4*Noscil_plot:end),'color',colours{1}); xlabel({'$t$'},'interpreter','latex'); ylabel({'$(V_i)_{FOM}$'},'interpreter','latex'); xlim([0 500]);

    nexttile(2); 
    ax2 = plot(rom.time,rom.x(:,1:4*Noscil_plot:end-3),'color',colours{2}); ylabel({'$(X_i)_{ROM}$'},'interpreter','latex'); xlim([0 500]); set(gca,'xtick',[])
    nexttile(5);
    plot(rom.time,rom.x(:,2:4*Noscil_plot:end-2),'color',colours{2}); ylabel({'$(Y_i)_{ROM}$'},'interpreter','latex'); xlim([0 500]); set(gca,'xtick',[])
    nexttile(8);
    plot(rom.time,rom.x(:,3:4*Noscil_plot:end-1),'color',colours{2}); ylabel({'$(Z_i)_{ROM}$'},'interpreter','latex'); xlim([0 500]); set(gca,'xtick',[])
    nexttile(11);
    plot(rom.time,rom.x(:,4:4*Noscil_plot:end),'color',colours{2}); xlabel({'$t$'},'interpreter','latex'); ylabel({'$(V_i)_{ROM}$'},'interpreter','latex'); xlim([0 500]);

    nexttile(3); 
    ax3 = plot(prom.time,prom.x(:,1:4*Noscil_plot:end-3),'color',colours{3}); ylabel({'$(X_i)_{pROM}$'},'interpreter','latex'); xlim([0 500]); set(gca,'xtick',[])
    nexttile(6);
    plot(prom.time,prom.x(:,2:4*Noscil_plot:end-2),'color',colours{3}); ylabel({'$(Y_i)_{pROM}$'},'interpreter','latex'); xlim([0 500]); set(gca,'xtick',[])
    nexttile(9);
    plot(prom.time,prom.x(:,3:4*Noscil_plot:end-1),'color',colours{3}); ylabel({'$(Z_i)_{pROM}$'},'interpreter','latex'); xlim([0 500]); set(gca,'xtick',[])
    nexttile(12);
    plot(prom.time,prom.x(:,4:4*Noscil_plot:end),'color',colours{3}); xlabel({'$t$'},'interpreter','latex'); ylabel({'$(V_i)_{pROM}$'},'interpreter','latex'); xlim([0 500]);
    
    lg = legend([ax1(1), ax2(1), ax3(1)], ["FOM", "$ROM$","$pROM$"], 'interpreter', 'latex');
    lg.Layout.Tile = 'north';

    % Visualize synchrony variable F/E(V^2)
    FIG2 = figure();
    t = tiledlayout(3,1);
    t.TileSpacing = 'tight';
    
    sgtitle('\textbf{Synchrony in each model}', 'FontSize', 40, 'interpreter', 'latex');
    
    nexttile;
    ax1 = plot(full.time,full.synchrony,'color',colours{1}); hold on; 
    ax2 = plot(rom.time,rom.synchrony,'color',colours{2}); 
    ax3 = plot(prom.time,prom.synchrony,'color',colours{3}); hold off;
    xlabel('$t$','interpreter','latex'); ylabel('$Q$','interpreter','latex'); set(gca,'xtick',[])
    title('Synchrony variable');
    
    % Visualize average gene concentration X(t)
    nexttile;
    plot(full.time,full.avg_gene,'color',colours{1}); hold on;
    plot(rom.time,rom.avg_gene,'color',colours{2}); 
    plot(prom.time,prom.avg_gene,'color',colours{3}); hold off;
    xlabel('$t$','interpreter','latex'); ylabel('$X$','interpreter','latex'); set(gca,'xtick',[])
    title('Average gene concentration');
    
    % Visualize order parameter
    nexttile;
    plot(full.time(end-n_asympt(1):end),full.order_param,'color',colours{1}); hold on;
    plot(rom.time(end-n_asympt(2):end),rom.order_param,'color',colours{2}); 
    plot(prom.time(end-n_asympt(3):end),prom.order_param,'color',colours{3}); hold off;
    xlabel('$t$','interpreter','latex'); ylabel('$\gamma$','interpreter','latex'); ylim([0 1]);
    title('Order parameter');
    
    lg = legend([ax1, ax2, ax3], ["FOM", "$ROM$","$pROM$"], 'interpreter', 'latex');
    lg.Layout.Tile = 'north';
    
end

