
clear; close all; clc;

set(0,'defaultAxesFontSize',35); set(0,'defaultLineLineWidth',3);

set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'CMU Serif');
set(0,'DefaultAxesFontName', 'CMU Serif');

addpath('functions')
addpath(genpath('../syst'))
addpath(genpath('figures'))

%% Definition of the systems

output = 'kuramoto_order';

factors = read_factorSpace('');
legend_cell = legendise(factors);

%% Reduction of the system and analysis

neurons = 100; 

    fprintf(['\n ## ROM for ', num2str(neurons), ' neurons \n']);
    pod_deim(neurons)
    
    
%% Control problem

neurons = 100; 

    fprintf(['\n ## ROM for ', num2str(neurons), ' neurons \n']);
    pod_deim(neurons)
