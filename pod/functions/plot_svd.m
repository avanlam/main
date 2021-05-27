function [FIG1, FIG2, FIG3] = plot_svd(U, S, U_total, V, T, V_total)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
set(0,'defaultAxesFontSize',35);
N = size(U{1}, 1);

%% Angle between the two subspaces

    theta = zeros(N,1);
    for n = 1:N
        theta(n) = subspace(U_total{1}(:,1:n),V_total{1}(:,1:n));
    end
            
    FIG1 = figure();   
    colours = colour_pairs('green');
    
    plot(1:N, theta, 'Marker', 'x','MarkerSize', 12, 'LineWidth', 3, 'Color', colours{2})
    ylim([0 pi/2]); yticks(0:pi/6:pi/2); yticklabels(["$$0$$","$$\pi/6$$","$$\pi/3$$","$$\pi/2$$"]);
    xlabel('$r/4$','interpreter','latex', 'FontSize', 40); 
    ylabel('$\theta \; (\mathrm{rad})$','interpreter','latex', 'FontSize', 40);
    title('\textbf{Angle between the ROM and pROM subspaces}', 'FontSize', 40, 'interpreter', 'latex');

    
%     r = 1;
%     t = linspace(0,2*pi);
%     
%     x = r*cos(t);
%     y = r*sin(t);
%     plot(x, y, ':', 'Color', colours{2}); hold on;
%     
%     t = linspace(0,theta);
%     x = r*cos(t);
%     y = r*sin(t);
%     fill([0,x,0],[0,y,0],colours{1})
%     plot([0,x,0],[0,y,0],'Color', colours{2})
%     text(1.1, 0.2, [num2str(theta),' rad'], 'FontSize', 30);
%     
%     xlim([-1.2 1.2]);ylim([-1.2 1.2]);
%     xticks([]);yticks([]);
%     title('\textbf{Angle between the ROM and pROM subspaces}', 'FontSize', 40, 'interpreter', 'latex');
% 
%     axis equal; axis off ;
    
%% Decay of the SVD eigenvalues

    FIG2 = figure();
    colours = colour_pairs('duo');

    t = tiledlayout(2,2);
    t.TileSpacing = 'tight';    
    
    titles = ["$X(t)$","$Y(t)$","$Z(t)$","$V(t)$"];
    for i = 1:4
        nexttile;
        semilogy(1:length(S{i}),S{i}/S{1}(1),'color',colours{2}, 'LineWidth', 3); hold on;
        semilogy(1:length(T{i}),T{i}/T{1}(1),'color',colours{1}, 'LineWidth', 3); hold off;
        title(titles(i));
        if i == 2 || i==4
            set(gca,'ytick',[])
        end
        if i == 1 || i==2
            set(gca,'xtick',[])
        end
        ylim([1e-20, 1]);
    end

    
    lg = legend(["$ROM$","$pROM$"], 'interpreter', 'latex');
    lg.Layout.Tile = 'north';
    xlabel(t, '$i$','interpreter','latex', 'FontSize', 30); 
    ylabel(t, '$\sigma_i/\sigma_1$','interpreter','latex', 'FontSize', 30);
    sgtitle('\textbf{Decay of the SVD}', 'FontSize', 40, 'interpreter', 'latex');

%% Comparison of the five first columns (most important vectors)
    
    FIG3 = figure();
    colours = colour_pairs('duo');

    t = tiledlayout(3,3);
    t.TileSpacing = 'tight';    
    
    for i = 1:9
        nexttile;
        
        plot(1:N,abs(U{1}(:,i)),'Marker', 'x', 'color',colours{2}); hold on;
        plot(1:N,abs(V{1}(:,i)),'Marker', 'x','color',colours{1}); hold off;
        ylim([0 0.2]); % yticks([0.099, 0.1 0.101]);
        title(['$\textbf{v}_', num2str(i), '$'],'interpreter','latex', 'FontSize', 30, 'Units', 'normalized', 'Position', [0.5, 0.8, 0]); %
        
        if i ~= 1 && i ~= 4 && i ~= 7
            set(gca,'ytick',[])
        end
        if i ~= 7 && i ~= 8 && i ~= 9
            set(gca,'xtick',[])
        end
        
    end
    
    lg = legend(["$ROM$","$pROM$"], 'interpreter', 'latex');
    lg.Layout.Tile = 'north';
    xlabel(t, '$i$','interpreter','latex', 'FontSize', 30); 
    sgtitle('\textbf{Comparison of the first POD modes}', 'FontSize', 40, 'interpreter', 'latex');
    
end

