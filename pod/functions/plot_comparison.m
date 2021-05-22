function FIG = plot_comparison(c_cell, name, method, outputs, legend_cell)

t_array = ["indices","ranking"];

%% Underline the most important according to Kuramoto

%[~, idx] = sort(sum(c_cell{1}, 2), 'descend');
t = c_cell{1}; [~, idx] = sort(t(:,3), 'descend');
clr_change = zeros(size(legend_cell));

legend_cell{end} = strcat('\', legend_cell{end});
for i = 1:5
    variab = legend_cell{idx(i)};
    legend_cell{idx(i)} = strcat('$\underline{', variab, '}$');
    clr_change(idx(i))=1;
end
for i = 6:16
    variab = legend_cell{idx(i)};
    legend_cell{idx(i)} = strcat('$', variab, '$');
end

%% Comparison figures
FIG = figure();
FIG.WindowState = 'fullscreen';

for i = 1:2
    subplot(2, 1, i);
    b = bar(c_cell{i});
    title( strcat("Sobol's ", t_array(i)) );
    xticks(1:length(legend_cell)); xticklabels(legend_cell); %xtickangle(90)

    % Pretty colours
    colours = [colour_pairs('blue'), colour_pairs('green'), colour_pairs('red')];
    b(1).FaceColor = colours{3};
    b(2).FaceColor = colours{1};
    b(3).FaceColor = colours{2};

    if i==1
        legend(outputs,'interpreter','latex', 'Location', 'best', 'FontSize', 20)
    else
        for j = 1:3
            xtips = b(j).XEndPoints;
            ytips = b(j).YEndPoints;
            labels = string(b(j).YData);
            text(xtips,ytips,labels,'HorizontalAlignment','center',...
                'VerticalAlignment','bottom', 'FontSize', 20)
        end
        yline(5,'b--')
    end
end

sgtitle(strcat("\textbf{Sobol's ", name, " according to different model indicators}"), 'FontSize', 45);
saveas(FIG,[method, '/comparison/sobol_',name,'_outputs'],'epsc')

%% Group comparison figures

% Pretty colours
colours = [colour_pairs('green'), colour_pairs('orange')];
    
FIG = figure();
FIG.WindowState = 'fullscreen';
x = 1:length(legend_cell); dx = 0.1;
dy = 0.1;

for i = 1:2
    y = median(c_cell{i},2);
    
    subplot(2,1,i)
    p1 = plot(repmat(x, 2,1), [min(c_cell{i}, [],2), max(c_cell{i}, [],2)]', 'Color', colours{1}, 'LineWidth', 4); hold on;
    p2 = plot(x, y, 'o', 'Color', colours{2}, 'MarkerSize', 10);
    if i==1
        y = round(y,2);
        legend([p1,p2], ["range", "median"],'interpreter','latex', 'Location', 'best', 'FontSize', 20);
    end
    text(x+dx, y, num2str(y), 'FontSize', 18)
    title( strcat("Sobol's ", t_array(i)) );
    xlim([0 length(legend_cell)]+0.5); xticks(1:length(legend_cell)); xticklabels(legend_cell);
    

end
sgtitle(strcat("\textbf{Sobol's ", name, " variation across the three model indicators}"), 'FontSize', 45);
saveas(FIG,[method,'/comparison/sobol_conv_',name,'_outputs'],'epsc')

%     b = boxplot(c_cell{i}', 'MedianStyle', 'target', 'Labels', legend_cell, 'Colors', [colours{2};colours{4}], 'Colorgroup', clr_change>0); hold on;
%     set(b,{'linew'},{2})
%     bp = gca; bp.YAxis.TickLabelInterpreter = 'latex'; bp.XAxis.TickLabelInterpreter = 'latex';

end