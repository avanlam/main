function [FIG] = plot_grid(param, name, x, y, lim, lb, ub)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
set(0,'defaultAxesFontSize',30);

% Write points in meshgrid format
[tmp1,tmp2]=meshgrid(x,y);
mesh_1 = reshape(tmp1',[],1);
mesh_2 = reshape(tmp2',[],1);
[XX,YY]=meshgrid(linspace(lb(1),ub(1),8*length(x)),linspace(lb(2),ub(2),8*length(y)));

    
FIG = figure();

t = tiledlayout(2,size(param,2));
t.TileSpacing = 'tight';

FIG.WindowState = 'fullscreen';
sgtitle(['\textbf{Kuramoto order parameter $$', name,'$$}'], 'interpreter', 'latex', 'FontSize', 45)        
titles = {'FOM', 'POD', 'pPOD', 'pDEIM'};
   
    for i = 1:size(param,2)
        nexttile;

        inter = scatteredInterpolant(mesh_1,mesh_2,param(:, i)); 
        surf(XX,YY,inter(XX, YY)); 
        colormap(parula); 
        view(0,90); shading interp; caxis(lim); 
        title(titles{i},'interpreter','latex');
        xlim([lb(1),ub(1)]); xticks(20:4:28); ylim([lb(2),ub(2)])
        axis square
        
        xlabel('$\Omega$','interpreter','latex', 'FontSize', 40); 
        if i==1
            ylabel('$L_0$','interpreter','latex', 'FontSize', 40);
%         else
%             set( gca, 'YTick', [] );
        end
        
    end
    
%     cb = colorbar('Location', 'north');
%     set(cb, 'Position', [.25 .83 .5 .0281])
    cb = colorbar; cb.Layout.Tile = 'north';
    cb.TickLabelInterpreter = 'latex';
    
    set(0,'defaultAxesFontSize',35);
end

