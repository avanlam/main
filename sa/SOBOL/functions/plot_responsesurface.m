
clear; close all; clc;

set(0,'defaultAxesFontSize',20); set(0,'defaultLineLineWidth',2);

set(0,'defaultAxesTickLabelInterpreter','latex');
set(0,'defaulttextinterpreter','latex');
set(0,'DefaultTextFontname', 'CMU Serif');
set(0,'DefaultAxesFontName', 'CMU Serif');

%% 

n_samples = 100;                            % 1. Size of the base sample: total number of function evaluation will be N * (numDim + 2)
seedNum = 123456789;                % 2. Seed for random number generator; put [] for autmatic randomization 
numN = 10;
funPath = '../../../syst';   % 4. Folder address: that includes model file, factor space
funFile = 'eval_circ';              % 5. Model/function file: MATLAB m-file without .m extension

outputs = {'sync_param','spectral_ampl','kuramoto_order'};

%% Generate the base sample from a unit hypercube
addpath(funPath);
factors = read_factorSpace(''); % read in the factor space

lb = factors.lb; ub = factors.ub; 
numDim = length(lb);  % number of factors

p = sobolset(numDim,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
p = scramble(p,'MatousekAffineOwen'); %MatousekAffineOwen scramling strategy. For no scrambling disable this line.
A = p(1:n_samples,:);

for j = 1 : n_samples
    A(j, :) = A(j, :) .* ( ub' - lb') + lb';
end
%% Run the function/model for all points in matrices A, B, and C and return the model response

yA = zeros(n_samples, 3);
for o = 1:length(outputs)
    output = outputs{o};
    for j = 1 : n_samples
        fprintf('Group run (base sample) #%g started. Running model %s %g times...\n', j, funFile, numDim + 2 );
        yA(j, o) = feval(funFile, [numN, A(j, :), o] );
    end
end

%% Plot the model response surface for all points in matrices A, B, and C
legend_cell = legendise(factors);
output_cell = {'$$\rho$$', '$$S$$', '$$R$$'};

c = [colour_groups('orange',5); colour_groups('green',5); colour_groups('blue',3); colour_groups('purple',2); colour_groups('red',3)];

for o = 1:length(outputs)
    output = outputs{o};
    
    FIG = figure();
    FIG.WindowState = 'fullscreen';
    sgtitle(['\textbf{Model response surface for} \verb|', output, '|'], 'interpreter','latex', 'FontSize', 40);
    for i = 1:numDim
        [x, idx] = sort(A(:,i));

        if i >15
            subplot(4,5, i+1)
        else
            subplot(4,5, i)
        end

        semilogy(x, yA(idx, o), '.', 'MarkerSize', 16, 'Color', c(i,:));
        xlim([lb(i), ub(i)]);
    %     yticks([1e-2, 1, 1e2, 1e4, 1e6]);
%         ylabel(output_cell(o), 'interpreter','latex', 'FontSize', 30);
        xlabel(legend_cell(i), 'interpreter','latex', 'FontSize', 30);
    end

    saveas(FIG, ['../figures/', output, '/response_surface'], 'epsc');
end