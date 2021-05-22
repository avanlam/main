function [FIG1,FIG2] = plot_sensitivity(funPath, ind1,ind2, labels)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

method = "Sobol Sensitivity";
title1 = "First Order";
title2 = "Total Order";

%% Bar Charts
FIG1 = figure();
FIG1.WindowState = 'fullscreen';
b = bar(1:numel(ind1), [abs(ind1)./sum(abs(ind1)); abs(ind2)]./sum(abs(ind2))); 
title("\textbf{Sensitivity of the parameters}", 'FontSize', 50);
xlabel("Parameter", 'FontSize', 40); ylabel(method, 'FontSize', 50);
ylim([0,1]); xticks(1:length(labels)); xticklabels(labels);
legend(title1,title2,'interpreter','latex', 'Location', 'best', 'FontSize', 40)

% Pretty colours
colours = colour_pairs('blue');
b(1).FaceColor = colours{1};
b(2).FaceColor = colours{2};

% Save figure full-screen
saveas(FIG1,strcat(funPath,'sobol_sens'),'epsc')

%%  Pie Charts
FIG2 = figure();
FIG2.WindowState = 'fullscreen';
sgtitle("\textbf{Importance of the parameters}", 'FontSize', 50);

ax1 = subplot(1,2,1);
pieHandle = pie(abs(ind1), labels);
title(title1)
set(findobj(ax1,'type','text'),'FontSize',40) 

pieAxis = get(pieHandle(1), 'Parent'); 
pieAxisPosition = get(pieAxis, 'Position');
newRadius = 1;   % Change the radius of the pie chart
newPieAxisPosition = pieAxisPosition .*[1 1 newRadius newRadius];
set(pieAxis, 'Position', newPieAxisPosition);

ax2 = subplot(1,2,2);
pieHandle = pie(abs(ind2), labels);
title(title2)
set(findobj(ax2,'type','text'),'FontSize',40) 

pieAxis = get(pieHandle(1), 'Parent'); 
pieAxisPosition = get(pieAxis, 'Position');
newRadius = 1;   % Change the radius of the pie chart
newPieAxisPosition = pieAxisPosition .*[1 1 newRadius newRadius];
set(pieAxis, 'Position', newPieAxisPosition);

% Save figure full-screen
set(FIG2, 'PaperPositionMode', 'auto')
set(FIG2, 'PaperSize', [41 29.7])
saveas(FIG2,strcat(funPath,'sobol_imp'),'pdf')

end

