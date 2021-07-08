% Inspired by: Saman Razavi in November 2019

clear; close all; clc;

set(0,'defaultAxesFontSize',35); 
set(0,'defaultAxesTickLabelInterpreter','latex');

set(0,'defaultLineLineWidth',2);

set(0,'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'CMU Serif');
set(0,'DefaultAxesFontName', 'CMU Serif');

addpath('functions')

%% Definition of the systems

addpath(genpath('figures'))

filename = 'VARS_inp.txt';
VARS_inp = read_VARS_inp(filename);

VARS_inp.neurons = 100; 
factors = read_factorSpace('');
VARS_inp.factors = factors;
legend_cell = legendise(factors);

%% Simulation of the systems

VARS_inp.outFldr = strcat(VARS_inp.outFldr, num2str(VARS_inp.neurons));

fprintf('A VARS experiment started: number of neurons = %g, number of stars = %g, minimum h = %g. \n', VARS_inp.neurons, VARS_inp.numStars, VARS_inp.grdSize);
VARS_out = VARS(VARS_inp);

figure(1); saveas(gcf,strcat(VARS_inp.outFldr,'/vars'),'fig')
figure(2); saveas(gcf,strcat(VARS_inp.outFldr,'/vars_conv'),'fig')
figure(3); saveas(gcf,strcat(VARS_inp.outFldr,'/vars_group'),'fig')

figure(1);
sgtitle('VARS comparison of Morris, Sobol and VARS', 'FontSize', 40);
xlabel('Parameters', 'FontSize', 30); ylabel('Parameter sensitivity', 'FontSize', 30);

GreenAsh = [160, 218, 169]./255
Mint = [0,161,112]./255
Marigold = [253,172,83]./255
