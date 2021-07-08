function [] = plot_presentation(U, U_total, S, V, V_total, T, n_asympt, full, rom, prom, param, x, y, lim, lb, ub, base)
%PLOT-PRESENTATION draws three different graphs explored in other
%functions, but adapted here in size, font and description for
%presentations

set(0,'defaultAxesFontSize',40);
N = size(U{1}, 1);

%% SVD

theta = zeros(N,1);
for n = 1:N
    theta(n) = subspace(U_total{1}(:,1:n),V_total{1}(:,1:n));
end

colours = colour_pairs('duo');
% colours{1} = [102,178,255]./255;
% colours{2} = [255,153,51]./255;
    
FIG = figure();

    t = tiledlayout(1,4);
    t.TileSpacing = 'tight';    
    
    nexttile;
    semilogy(1:length(S{1}),S{1}/S{1}(1),'color',colours{2}, 'LineWidth', 4); hold on;
    semilogy(1:length(T{1}),T{1}/T{1}(1),'color',colours{1}, 'LineWidth', 4); hold off;
    xlabel('$i$','interpreter','latex', 'FontSize', 30); 
    ylim([1e-20, 1]); title('$\mathbf{\lambda}_i$','interpreter','latex', 'FontSize', 45);
    legend('ROM', 'pROM', 'interpreter','latex', 'FontSize', 30);
    
    nexttile;
	plot(1:N,U{1}(:,2), 'color',colours{2}); hold on;
	plot(1:N,V{1}(:,2), ':','color',colours{1}); hold off;
    xlabel('$i$','interpreter','latex', 'FontSize', 30); 
    ylim([-0.2 0.2]); % yticks([0.099, 0.1 0.101]);
    title('$\textbf{v}_2$','interpreter','latex', 'FontSize', 45); %
        
    nexttile;
	plot(1:N,-U{1}(:,5), 'color',colours{2}); hold on;
	plot(1:N,V{1}(:,5), ':','color',colours{1}); hold off;
	xlabel('$i$','interpreter','latex', 'FontSize', 30); 
    ylim([-0.2 0.2]); % yticks([0.099, 0.1 0.101]);
	title('$\textbf{v}_5$','interpreter','latex', 'FontSize', 45); %
                
            
    nexttile;
    plot(1:20, theta(1:20), 'Marker', 'x', 'Color', colours{2}, 'LineWidth', 4, 'MarkerSize', 10)
    ylim([0 pi/2]); yticks(0:pi/6:pi/2); yticklabels(["$$0$$","$$\pi/6$$","$$\pi/3$$","$$\pi/2$$"]);
    xlabel('$r/4$','interpreter','latex', 'FontSize', 30); 
    title('$\mathbf{\theta}$', 'FontSize', 45, 'interpreter', 'latex');

pause;
exportgraphics(FIG,strcat('figures/presentation/svd_vectors.eps'),'ContentType','vector','BackgroundColor','none')


%% Basal value

if numel(n_asympt) == 1
    n_asympt = [n_asympt n_asympt n_asympt];
end

colours = colour_pairs('spring');

    % Visualize synchrony variable F/E(V^2)
    FIG2 = figure();
    t = tiledlayout(1,3);
    t.TileSpacing = 'tight';
    
    % Visualize oscillators
    nexttile; Noscil_plot = N;
    plot(full.time(end-n_asympt(1):end),full.x(end-n_asympt(1):end,1:4*Noscil_plot:end-3),'color',colours{1}); hold on;
    plot(rom.time(end-n_asympt(2):end),rom.x(end-n_asympt(1):end,1:4*Noscil_plot:end-3),'x','color',colours{2}, 'MarkerSize', 12); 
    plot(prom.time(end-n_asympt(3):end),prom.x(end-n_asympt(1):end,1:4*Noscil_plot:end-3),'.','color',colours{3}, 'MarkerSize', 16); hold off;
    ylim([0 0.23]); title('$\mathbf{X}_1$','interpreter','latex');
   
    % Visualize average gene concentration X(t)
    nexttile;
    plot(full.time(end-n_asympt(1):end),full.avg_gene(end-n_asympt(1):end),'color',colours{1}); hold on;
    ax2 = plot(rom.time(end-n_asympt(2):end),rom.avg_gene(end-n_asympt(2):end),'x','color',colours{2}, 'MarkerSize', 12); 
    plot(prom.time(end-n_asympt(3):end),prom.avg_gene(end-n_asympt(3):end),'.','color',colours{3}, 'MarkerSize', 16); hold off;
    xlabel('$t$','interpreter','latex');
    ylim([0 0.23]); title('$\mathbf{\bar X}$','interpreter','latex');
    
    % Visualize order parameter
    nexttile;
    ax1 = plot(full.time(end-n_asympt(1):end),full.order_param,'color',colours{1}); hold on;
    plot(rom.time(end-n_asympt(2):end),rom.order_param,'.','color',colours{2}, 'MarkerSize', 12); 
    ax3 = plot(prom.time(end-n_asympt(3):end),prom.order_param,'.','color',colours{3}, 'MarkerSize', 16); hold off; 
    ylim([0.5 1]); title('$\mathbf{\gamma}$','interpreter','latex');
    legend([ax1, ax2, ax3], "FOM", "$ROM$","$pROM$", 'interpreter', 'latex', 'location', 'best');

    pause;
    exportgraphics(FIG2,strcat('figures/presentation/basal.eps'),'ContentType','vector','BackgroundColor','none')
    
%% PARAMETRIC VARIATION
    
% Write points in meshgrid format
[tmp1,tmp2]=meshgrid(x,y);
mesh_1 = reshape(tmp1',[],1);
mesh_2 = reshape(tmp2',[],1);
[XX,YY]=meshgrid(linspace(lb(1),ub(1),8*length(x)),linspace(lb(2),ub(2),8*length(y)));
base_sample_rom = 0.5.* ( ub' - lb') + lb';
    
FIG3 = figure();

t = tiledlayout(2,3);
t.TileSpacing = 'none';
title(t, '$L_2$ \textbf{error}','interpreter','latex', 'FontSize', 45);

titles = {'ROM', 'pROM'};
   
for i = 1:size(param,2)
    nexttile;

    inter = scatteredInterpolant(mesh_1,mesh_2,param(:, i)); 
    surf(XX,YY,inter(XX, YY)); 
    colormap(parula); 
    view(0,90); shading interp; caxis(lim); 
    title(titles{i},'FontSize', 40, 'interpreter','latex');
    xlim([lb(1),ub(1)]); xticks(20:4:28); ylim([lb(2),ub(2)])
    axis square

    if i==1
        ylabel('$L_0$','interpreter','latex', 'FontSize', 40);
    elseif i ==2
        xlabel('$\Omega$','interpreter','latex', 'FontSize', 40); 
    end
end
    
    cb = colorbar; cb.Layout.Tile = 'north';
    cb.TickLabelInterpreter = 'latex';
    
    nexttile;
    s(2) = scatter(base(:,1),base(:,2), 300, colours{3}, 'filled', 'Marker', 's'); hold on;
    s(1) = scatter(base_sample_rom(:,1),base_sample_rom(:,2), 300, colours{2}, 'filled', 'Marker', 's'); hold off
    xlim([lb(1),ub(1)]); ylim([lb(2),ub(2)]); axis square
    
    [~, objh] = legend(s,'ROM','pROM', 'interpreter','latex', 'FontSize', 25, 'Location', 'northeast');
    objh(3).Children.MarkerSize = 16; objh(4).Children.MarkerSize = 16;    

    pause;
    exportgraphics(FIG3,strcat('figures/presentation/parametric.eps'),'ContentType','vector','BackgroundColor','none')
    
%% revert
  
set(0,'defaultAxesFontSize',35);

end

