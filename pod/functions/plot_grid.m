function [FIG] = plot_grid(output, name, x, y, lim, lb, ub)
%PLOT_GRID draws the synchronisation measures of the ROM, pROM and FOM when
% two constitutive parameters are varied simultaneously across a grid
% INPUT :
%   * 'output' is a matrix containing the output to be drawn, where each
%   column corresponds to a different model (ROM, pROM and FOM)
%   * 'name' contains the name of the synchronisation output plotted
%   * 'x' and 'y' gives the coordinates of the parameter grid
%   * 'lb' and 'ub' indicate the lower and upper bounds on the parameter
%   grid, which will define the ends of the graph
%   * 'lim' indicates the bounds on the output's value, which will define
%   the ends of the colorbar

set(0,'defaultAxesFontSize',30);

% Write points in meshgrid format
[tmp1,tmp2]=meshgrid(x,y);
mesh_1 = reshape(tmp1',[],1);
mesh_2 = reshape(tmp2',[],1);
[XX,YY]=meshgrid(linspace(lb(1),ub(1),8*length(x)),linspace(lb(2),ub(2),8*length(y)));

    
FIG = figure();

t = tiledlayout(2,size(output,2));
t.TileSpacing = 'tight';

FIG.WindowState = 'fullscreen';
sgtitle(['\textbf{Kuramoto order parameter $$', name,'$$}'], 'interpreter', 'latex', 'FontSize', 45)        
titles = {'FOM', 'POD', 'pPOD', 'pDEIM'};
   
    for i = 1:size(output,2)
        nexttile;

        inter = scatteredInterpolant(mesh_1,mesh_2,output(:, i)); 
        surf(XX,YY,inter(XX, YY)); 
        colormap(parula); 
        view(0,90); shading interp; caxis(lim); 
        title(titles{i},'interpreter','latex');
        xlim([lb(1),ub(1)]); xticks(20:4:28); ylim([lb(2),ub(2)])
        axis square
        
        xlabel('$\Omega$','interpreter','latex', 'FontSize', 40); 
        if i==1
            ylabel('$L_0$','interpreter','latex', 'FontSize', 40);
        end
        
    end
    
    cb = colorbar; cb.Layout.Tile = 'north';
    cb.TickLabelInterpreter = 'latex';
    
    set(0,'defaultAxesFontSize',35);
end

