% ****         Variogram Analysis of Response Surfaces (VARS)         ****
% ****                              for                               ****
% ****           Global Sensitivity and Uncertainty Analysis          ****
% ****                      © Saman Razavi 2018                       ****
% ****                                                                ****
% ****        Version 1 written by Saman Razavi in 2013-2015          ****
% ****                                                                ****
% ****       Updated to Version 2 by Saman Razavi in 2017-2018        ****
% **** enabled with dynamic SA, grouping, plotting, etc. capabilities ****
% ****                                                                ****
% ****    Updated to Version 2.1 by Saman Razavi in November 2019     ****
% ****    with improved plotting and some other small improvements    ****
% ****      \                                                 /       ****
% ****       --- THIS IS THE MAIN FUNCTION TO EXECUTE VARS ---        ****
% ************************************************************************
clear; close all; clc;

set(0,'defaultAxesFontSize',35); set(0,'defaultLineLineWidth',2);

set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'CMU Serif');
set(0,'DefaultAxesFontName', 'CMU Serif');

addpath('functions')

%% Definition of the systems

method = 'figures_determ'; addpath(genpath(method))

filename = 'VARS_inp.txt';
VARS_inp = read_VARS_inp(filename);

VARS_inp.neurons = 100; 
VARS_inp.outputs = 1; output_names = 'kuramoto_order';

factors = read_factorSpace('');
VARS_inp.factors = factors;
legend_cell = legendise(factors);

%% Simulation of the systems

VARS_inp.outFldr = strcat(VARS_inp.outFldr, num2str(VARS_inp.neurons));

VARS_inp.factors
fprintf('A VARS experiment started: number of neurons = %g, number of stars = %g, minimum h = %g. \n', VARS_inp.neurons, VARS_inp.numStars, VARS_inp.grdSize);
VARS_out = VARS(VARS_inp);

figure(1); saveas(gcf,strcat(VARS_inp.outFldr,'/vars_', num2str(VARS_inp.neurons)),'fig')
figure(2); saveas(gcf,strcat(VARS_inp.outFldr,'/vars_conv_', num2str(VARS_inp.neurons)),'fig')
figure(3); saveas(gcf,strcat(VARS_inp.outFldr,'/vars_group_', num2str(VARS_inp.neurons)),'fig')


% %% Simulation of the systems
% for o = 1:length(outputs)
%     output = outputs(o);
%     for n = 1:length(neurons)
%         num_n = neurons(n);
% 
%         VARS_copy = VARS_inp;
%         
%         VARS_copy.num_n = num_n;
%         VARS_copy.outFldr = strcat(VARS_copy.outFldr, num2str(num_n));
% 
%         VARS_copy.factors.ub = [num_n, VARS_copy.factors.ub', output]';
%         VARS_copy.factors.lb = [num_n, VARS_copy.factors.lb', output]';
%         VARS_copy.factors.name = [{'num_n'}, VARS_copy.factors.name', {'output'}]';
%         VARS_copy.factors.numDim = 18;
%         VARS_copy.numDim = 18;
%                 
%         VARS_copy.factors
%         fprintf('A VARS experiment started: number of neurons = %g, number of stars = %g, minimum h = %g. \n', VARS_copy.num_n, VARS_inp.numStars, VARS_inp.grdSize);
%         VARS_out = VARS(VARS_copy);
%         close all;
% 
%     end
% end

%% Document the evolution
results_all = struct();

indices = zeros(length(neurons), length(legend_cell), 4);
lb = zeros(length(neurons), length(legend_cell), 2);
ub = zeros(length(neurons), length(legend_cell), 2);

for o = 1:length(outputs)
    output = outputs{o};
    for n = 1:length(neurons)
        numN = neurons(n);

        load([output,'/Results_Sobol_', num2str(numN)]);
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

