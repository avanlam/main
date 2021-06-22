function [FIG] = plot_grid_error(param, name, x, y, lim, lb, ub, base)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% Write points in meshgrid format
[tmp1,tmp2]=meshgrid(x,y);
mesh_1 = reshape(tmp1',[],1);
mesh_2 = reshape(tmp2',[],1);
[XX,YY]=meshgrid(linspace(lb(1),ub(1),8*length(x)),linspace(lb(2),ub(2),8*length(y)));
base_sample_rom = 0.5.* ( ub' - lb') + lb';
    
FIG = figure();

t = tiledlayout(2,3);
t.TileSpacing = 'tight';

FIG.WindowState = 'fullscreen';
sgtitle(['\textbf{Parametric variation on the $L_2$ error}'], 'interpreter', 'latex', 'FontSize', 50)        
titles = {'ROM', 'pROM'};
   
    for i = 1:size(param,2)
        nexttile;

        inter = scatteredInterpolant(mesh_1,mesh_2,param(:, i)); 
        surf(XX,YY,inter(XX, YY)); 
        colormap(parula); 
        view(0,90); shading interp; caxis(lim); 
        title(['$', name, ' $(', titles{i}, ')'],'interpreter','latex');
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
    colours = colour_pairs('spring');
    s(2) = scatter(base(:,1),base(:,2), 300, colours{3}, 'filled', 'Marker', 's'); hold on;
    s(1) = scatter(base_sample_rom(:,1),base_sample_rom(:,2), 300, colours{2}, 'filled', 'Marker', 's'); hold off
    xlim([lb(1),ub(1)]); ylim([lb(2),ub(2)]); axis square
    
    [~, objh] = legend(s,'ROM','pROM', 'interpreter','latex', 'FontSize', 25, 'Location', 'northeast');
    objh(3).Children.MarkerSize = 16; objh(4).Children.MarkerSize = 16;    

end

