function [ Sobol_inp , Sobol_out ] = sobol(method, numN, name, output, color)
% ****       Sobol' Variance-based Global Sensitivity Analysis        ****
% ****     This numerical (Monte Carlo-based) integration method      ****
% ****            is adopted from  Saltelli et al. (2008)             ****
% ****  where page 164 states "This procedure is the best available   ****
% **** today for computing indices based purely on model evaluations" ****
% ****                 Written by © Saman Razavi 2018                 ****
% ************************************************************************
%% Control Parameters and Specifications
N = 500;                                   % 1. Size of the base sample: total number of function evaluation will be N * (numDim + 2)
seedNum = 123456789;                       % 2. Seed for random number generator; put [] for autmatic randomization 
funPath = ['../syst_', method];            % 3. Folder address: that includes model file, factor space
funFile = ['eval_circ_', name];            % 4. Model/function file: MATLAB m-file without .m extension
smplMtd = 'LHS';                           % 5. Sampling Method: RND, LHS, SymLHS, PLHS, SobolSeq, or Halton for generation of star centers; if blank, default is LHS
bootstrapFlag = 1;                         % 6. Bootstrap Flag: enter 1 to bootstrap, 0 not to bootstrap
bootstrapSize = 1000;                      % 7. Bootstrap Size: number of sampling iterations with replacement - if bootstrap flag = 0, this will be ignored.
confdLvl = 0.9;                            % 8. Confidence Level: to report bootstrap-based lower and upper bounds on VARS results; if bootstrapFlag = 0, this line will be ignored.
numGrp = 3;                                % 9. User-specified number of groups: if blank, VARS-TOOL will find the optimal number; if bootstrapFlag = 0, this line will be ignored.
plotFlag = 0;                              % 10. Plot Flag: enter 1 to plot, 0 not to plot
chosen_output = find(["sync_full","spectral_ampl","kuramoto_order"]==output);
%% Store Sobol Control Parameters and Specifications
Sobol_inp.N = N; 
Sobol_inp.seedNum = seedNum;
Sobol_inp.funFile = funFile;
Sobol_inp.funPath = funPath;
Sobol_inp.funFile = funFile;
Sobol_inp.smplMtd = smplMtd;
%% Randomization
if isempty( seedNum ) == false
    rand('state',seedNum); 
end
%% Generate the base sample from a unit hypercube
addpath(funPath);
factors = read_factorSpace(funPath, name); % read in the factor space
lb = factors.lb; ub = factors.ub; 
numDim = length(lb);  % number of factors
switch smplMtd
    case 'RND'
        baseSample = rand(N, numDim * 2); %uniformly distributed random numbers
    case 'LHS'
        baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
    case 'SymLHS'
        baseSample = SymLHS( N , numDim * 2,'maxmin', 10 );
    case 'PLHS'
        numSlices = ceil(N/sliceSize);
        if numSlices > 1
            smpSize = numSlices * sliceSize;
            p = PLHSdesign(smpSize, numDim * 2, numSlices, 10, 'maxmin');
            baseSample = p(1:N,:);
            % generates optimal PLHS based on 'maxmin' criterion.
            % by default the number of iteration is 10.
        else
            baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
        end
    case 'Sobolseq'
        p = sobolset(numDim * 2,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
        p = scramble(p,'MatousekAffineOwen'); %MatousekAffineOwen scramling strategy. For no scrambling disable this line.
        baseSample = p(1:N,:);
    case 'Halton'
        p = haltonset(numDim * 2,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
        p = scramble(p,'RR2'); %PR2 scramling strategy. For no scrambling disable this line.
        baseSample = p(1:N,:);
    otherwise
        baseSample = lhsdesign(N, numDim * 2,'criterion','maximin', 'iterations',100 );
end
%% Define sub-sample matrices A, B, and C, and scale them onto the actual factor space
A = baseSample(:, 1 : numDim);
B = baseSample(:, numDim + 1 : end);
for j = 1 : N
    A(j, :) = A(j, :) .* ( ub' - lb') + lb';
    B(j, :) = B(j, :) .* ( ub' - lb') + lb';
end
for i = 1 : numDim
    C{i} = B;
    C{i}(:, i) = A(:, i);
end

%% Run the function/model for all points in matrices A, B, and C and return the function/model response
[yA, yB, yC ] = func_eval(A, B, C, funPath, funFile, numN, chosen_output);
%% Calculate Sobol's First-Order (FO) and Total-Order (TO) sensitivity indices and rank factors based on TO
[ FO, TO, V, Mu ] = sobol_calc (yA, yB, yC);
rnkFO = factor_ranking(FO);
rnkTO = factor_ranking(TO);
Sobol_out.FO = FO;
Sobol_out.TO = TO;
Sobol_out.rnkFO = rnkFO;
Sobol_out.rnkTO = rnkTO;
Sobol_out.V = V;
Sobol_out.Mu = Mu;

%% Bootstrap
if bootstrapFlag == 1
    [ Sobol_out.bootstrap ] = bootstrap_sobol (yA, yB, yC, bootstrapSize, confdLvl, rnkFO, rnkTO, numGrp);
end
%% Store results
save ([method, '/', output, '/Results_Sobol_', name, '_', num2str(numN)], 'Sobol_inp', 'Sobol_out');
%% Plot results

if plotFlag == 1

    set(0,'defaultAxesFontSize',35); set(0,'defaultLineLineWidth',2);

    set(0,'defaultAxesTickLabelInterpreter','latex');
    set(0,'defaulttextinterpreter','latex');
    set(0,'DefaultTextFontname', 'CMU Serif');
    set(0,'DefaultAxesFontName', 'CMU Serif');

    legendCell = factors.name;
    for i = 1:numel(legendCell)
        legendCell{i} = strcat('$$', legendCell{i},'$$');
    end
    legendCell{end} = "$$\Omega$$";

    plot_sensitivity(strcat(funPath, "/Figures/"), strcat(name, '_' ,num2str(numN)), Sobol_out.FO, Sobol_out.TO, legendCell, color);
end

end