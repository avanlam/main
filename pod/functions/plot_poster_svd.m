function [] = plot_poster_svd(U, V, N)
%PLOT-POSTER-SVD draws the singular valeu decomposition in a graph adapted
% in terms of size and font to the scientific poster

colours = colour_pairs('duo');
colours{1} = [102,178,255]./255;
colours{2} = [255,153,51]./255;
set(0,'defaultAxesFontSize',45);
    
FIG = figure();

    t = tiledlayout(1,2);
    t.TileSpacing = 'tight';    
    
    nexttile;
               
	plot(1:N,-U{1}(:,2),'Marker', 'x', 'color',colours{2}); hold on;
	plot(1:N,V{1}(:,2),'Marker', 'x','color',colours{1}); hold off;
    ylim([-0.2 0.2]); % yticks([0.099, 0.1 0.101]);
    title('$\textbf{v}_2$','interpreter','latex', 'FontSize', 60); %
        
    nexttile;
                
	plot(1:N,-U{1}(:,5),'Marker', 'x', 'color',colours{2}); hold on;
	plot(1:N,V{1}(:,5),'Marker', 'x','color',colours{1}); hold off;
	xlabel(t, '$i$','interpreter','latex', 'FontSize', 30); 
    ylim([-0.2 0.2]); % yticks([0.099, 0.1 0.101]);
	title('$\textbf{v}_5$','interpreter','latex', 'FontSize', 60); %
                
pause;
exportgraphics(FIG,strcat('figures/poster_svd.eps'),'ContentType','vector','BackgroundColor','none')

set(0,'defaultAxesFontSize',35);
end