
clear; close all; clc;

set(0,'defaultAxesFontSize',35); set(0,'defaultLineLineWidth',2);

set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'CMU Serif');
set(0,'DefaultAxesFontName', 'CMU Serif');

addpath('functions')
addpath(genpath('figures'))
addpath(genpath('../../syst'))

%% Definition of the systems

neurons = [3, 5, 10, 20, 50, 100, 200, 500, 1000];
outputs = {'sync_param','spectral_ampl','kuramoto_order'};
names = {'$\rho$', '$S$', '$\gamma$'};
simulFlag = 0;

factors = read_factorSpace('');
legend_cell = legendise(factors);

%% Simulation of the systems

if simulFlag
    for o = 1:length(outputs)
        output = outputs{o};
        for n = 1:length(neurons)
            numN = neurons(n);

                sobol(numN, output);
                close all;

        end
    end
end

%% Document the evolution
results_all = struct();

indices = zeros(length(neurons), length(legend_cell), 4);
lb = zeros(length(neurons), length(legend_cell), 2);
ub = zeros(length(neurons), length(legend_cell), 2);

for o = 1:length(outputs)
    output = outputs{o};
    for n = 1:length(neurons)
        numN = neurons(n);

        load([output,'/SobolSeq18/Results_Sobol_', num2str(numN)]);
        results_all.(['r_', num2str(numN)]) = Sobol_out;
   
        indices(n, :, 1, o) = results_all.(['r_', num2str(numN)]).FO;
        indices(n, :, 2, o) = results_all.(['r_', num2str(numN)]).TO;
        indices(n, :, 3, o) = results_all.(['r_', num2str(numN)]).rnkFO;
        indices(n, :, 4, o) = results_all.(['r_', num2str(numN)]).rnkTO;
        
        lb(n, :, 1, o) = results_all.(['r_', num2str(numN)]).bootstrap.FO_low;
        lb(n, :, 2, o) = results_all.(['r_', num2str(numN)]).bootstrap.TO_low;
        
        ub(n, :, 1, o) = results_all.(['r_', num2str(numN)]).bootstrap.FO_upp;
        ub(n, :, 2, o) = results_all.(['r_', num2str(numN)]).bootstrap.TO_upp;
    end

end

%% Evolution of Sobol sensitivity analysis
fig_cell = {'sobol_fo_evolution', 'sobol_to_evolution', 'sobol_fo_rnk_evolution', 'sobol_to_rnk_evolution'};
tit_cell = {'FO indices', 'TO indices', 'FO ranking', 'TO ranking'};

for o = 1:length(outputs)
    
    plot_evolution_together('sobol_fo', cat(3,indices(:,:,1,o), lb(:,:,1,o), ub(:,:,1,o)), indices(:,:,3,o), neurons, legend_cell, {'FO indices', 'FO ranking'}, outputs{o}, names{o})
    plot_evolution_together('sobol_to', cat(3,indices(:,:,2,o), lb(:,:,2,o), ub(:,:,2,o)), indices(:,:,4,o), neurons, legend_cell, {'TO indices', 'TO ranking'}, outputs{o}, names{o})
    close all
end

% plot_evolution_poster(cat(3,indices(:,:,1,o), lb(:,:,1,o), ub(:,:,1,o)), neurons, legend_cell, outputs{3})
% plot_presentation(cat(3,indices(:,:,1,o), lb(:,:,1,o), ub(:,:,1,o)), neurons, legend_cell, {ind_fo, rnk_fo}, output_array)

%% Comparison between outputs

ind_fo = reshape(indices(end, :, 1, :), [length(legend_cell), length(outputs)]);
ind_to = reshape(indices(end, :, 2, :), [length(legend_cell), length(outputs)]);
rnk_fo = reshape(indices(end, :, 3, :), [length(legend_cell), length(outputs)]);
rnk_to = reshape(indices(end, :, 4, :), [length(legend_cell), length(outputs)]);

output_array = ["\verb|sync_param| ($\rho$)", "\verb|spectral_ampl| ($S$)", "\verb|kuramoto_order| ($\gamma$)"];
plot_comparison({ind_fo, rnk_fo}, 'FO', output_array, factors.name);
plot_comparison({ind_to, rnk_to}, 'TO', output_array, factors.name);

close all

%% 

% sobol(500, 'kuramoto_order');
load('kuramoto_order/presentation/Results_Sobol_500');
load('kuramoto_order/presentation/basesample');

colours = colour_pairs('spring');
FIG = figure;
    ax = gca();
    h = pie(ax, Sobol_out.FO.*2, legend_cell);
    ax.Colormap = reshape(cell2mat(colours), 3, 3)'; 

    h(2).FontSize = 65;
    h(4).FontSize = 65;
    h(6).FontSize = 65;
    
exportgraphics(FIG,strcat('figures/presentation/sobol_pie.eps'),'ContentType','vector','BackgroundColor','none')

FIG = figure();
   FIG.WindowState = 'fullscreen'; 
   scatter3(A(1:50,1), A(1:50,2), A(1:50,3), 300, colours{1}, 'filled', 'Marker', 's'); hold on
   xlim([lb(1),ub(1)]); ylim([lb(2),ub(2)]); axis square
   xlabel('$K_0$','interpreter','latex', 'FontSize', 40); 
   ylabel('$L_0$','interpreter','latex', 'FontSize', 40);
   zlabel('$\Omega$','interpreter','latex', 'FontSize', 40);

exportgraphics(FIG,strcat('figures/presentation/sobol_input.eps'),'ContentType','vector','BackgroundColor','none')

   
   for i = 1:3
       FIG = figure();
       scatter(1:50, A(1:50,i), 300, colours{i}, 'filled', 'Marker', 's'); hold on
       ylim([factors.lb(i), factors.ub(i)]); axis square; set(gca, 'XTick', []); set(gca, 'YTick', []);
       title(legend_cell{i},'interpreter','latex', 'FontSize', 65);
       exportgraphics(FIG,strcat('figures/presentation/sobol_input_', num2str(i),'.eps'),'ContentType','vector','BackgroundColor','none')
   end

