
clear; close all; clc;

set(0,'defaultAxesFontSize',16); set(0,'defaultLineLineWidth',2);

set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'CMU Serif');
set(0,'DefaultAxesFontName', 'CMU Serif');

%% 

N = 100;                            % 1. Size of the base sample: total number of function evaluation will be N * (numDim + 2)
seedNum = 123456789;                % 2. Seed for random number generator; put [] for autmatic randomization 
numN = 10;
funPath = '../../multiN';              % 4. Folder address: that includes model file, factor space
funFile = 'eval_circ_all';          % 5. Model/function file: MATLAB m-file without .m extension
name = 'all';
output = 'kuramoto_order';
chosen_output = find(["sync_full","spectral_ampl","kuramoto_order"]==output);

%% Generate the base sample from a unit hypercube
addpath(funPath);
factors = read_factorSpace(funPath, name); % read in the factor space
lb = factors.lb; ub = factors.ub; 
numDim = length(lb);  % number of factors
baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
A = baseSample(:, 1 : numDim);

%% Run the function/model for all points in matrices A, B, and C and return the model response
currDir = pwd;
cd(funPath);

yA = zeros(N, 1);
for j = 1 : N
    fprintf('Group run (base sample) #%g started. Running model %s %g times...\n', j, funFile, numDim + 2 );
    yA(j, 1) = feval(funFile, [numN, A(j, :), chosen_output] );
end
cd (currDir);

%% Plot the model response surface for all points in matrices A, B, and C
legendCell = factors.name;
for i = 1:numel(legendCell)
    legendCell{i} = strcat('$$', legendCell{i},'$$');
end
legendCell{end} = "$$\Omega$$";
c = [colour_groups('green',4); colour_groups('blue',4); colour_groups('purple',4); colour_groups('red',4)];

FIG = figure();
FIG.WindowState = 'fullscreen';
sgtitle(['Model response surface for \verb|', output, '|'], 'interpreter','latex', 'FontSize', 30);
for i = 1:numDim
    [x, idx] = sort(A(:,i));
    
    subplot(4,4, i)
    semilogy(x,yA(idx), 'Color', c(i,:));
    %yticks([1e-2, 1, 1e2, 1e4, 1e6]);
    title(legendCell(i), 'interpreter','latex', 'FontSize', 20,'fontweight','bold');
end