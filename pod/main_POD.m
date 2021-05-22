
clear; close all; clc;

set(0,'defaultAxesFontSize',35); set(0,'defaultLineLineWidth',2);

set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'CMU Serif');
set(0,'DefaultAxesFontName', 'CMU Serif');

addpath('functions')

%% Definition of the systems

method = 'determ';
addpath(genpath(['../syst_', method]))
addpath(genpath(['figures_', method]))

output = 'kuramoto_order'; name = 'all'; color = 'green';

factors = read_factorSpace('');

legendCell = factors.name;
for i = 1:numel(legendCell)
    legendCell{i} = strcat('$$', legendCell{i},'$$');
end
legendCell{end} = "$$\Omega$$";

%% Reduction of the system and analysis

neurons = 100; 

    fprintf(['\n ## ROM for ', num2str(neurons), ' neurons \n']);
    pod_deim(method, neurons)


