function [ Morris_inp , Morris_out ] = main_Morris()
% ****      Morris Derivative-based Global Sensitivity Analysis       ****
% ****   This numerical implementation computes elementary-effects    ****
% ****      indices from Morris (1991), Campolongo et al. (2007)      ****
% ****                and  Sobol and Kucherenko (2009)                ****
% ****                                                                ****
% ****                 Written by © Saman Razavi 2018                 ****
% ************************************************************************
%% Control Parameters and Specifications
numStrt = 100;                        % 1. Size of the base sample: total number of function evaluation will be numStrt * (numDim + 1)
stpSize = 0.1;                        % 2. Step size: for computing numerical derivatives; factor ranges will be scaled to [0-1]
seedNum = 123456789;                  % 3. Seed for random number generator; put [] for autmatic randomization
funPath = 'C:\2019_11_16_VARS-Tool-v2.1\HBV-SASK'; % 4. Folder address: that includes model file, factor space
funFile = 'eval_HBV_SASK';            % 5. Model/function file: MATLAB m-file without .m extension
smplMtd = 'LHS';                      % 6. Sampling Method: RND, LHS, SymLHS, PLHS, SobolSeq, or Halton for generation of star centers; if blank, default is LHS
sliceSize = 10;                       % 7. Slice Size: needed only for PLHS sampling strategy
bootstrapFlag = 1;                    % 8. Bootstrap Flag: enter 1 to bootstrap, 0 not to bootstrap
bootstrapSize = 1000;                 % 9. Bootstrap Size: number of sampling iterations with replacement - if bootstrap flag = 0, this will be ignored.
confdLvl = 0.9;                       % 10. Confidence Level: to report bootstrap-based lower and upper bounds on VARS results; if bootstrapFlag = 0, this line will be ignored.
numGrp = 3;                           % 11. User-specified number of groups: if blank, VARS-TOOL will find the optimal number; if bootstrapFlag = 0, this line will be ignored.
%% Store Sobol Control Parameters and Specifications
Morris_inp.numStrt = numStrt;
Morris_inp.seedNum = seedNum;
Morris_inp.funFile = funFile;
Morris_inp.funPath = funPath;
Morris_inp.funFile = funFile;
Morris_inp.smplMtd = smplMtd;
%% Randomization
if isempty( seedNum ) == false
    rand('state',seedNum);
end
%% Generate the start points of chains from a unit hypercube
factors = read_factorSpace(funPath); % read in the factor space
lb = factors.lb; ub = factors.ub;
numDim = length(lb);  % number of factors
MorrisCost = numStrt * (numDim + 1); % total number of function evaluations (model runs)
switch smplMtd
    case 'RND'
        strtPoints = rand(numStrt, numDim ); %uniformly distributed random numbers
    case 'LHS'
        strtPoints = lhsdesign(numStrt, numDim, 'criterion', 'maximin', 'iterations', 100 );
    case 'SymLHS'
        strtPoints = SymLHS( numStrt , numDim, 'maxmin', 10 );
    case 'PLHS'
        numSlices = ceil(numStrt/sliceSize);
        if numSlices > 1
            smpSize = numSlices * sliceSize;
            p = PLHSdesign(smpSize, numDim, numSlices, 10, 'maxmin');
            strtPoints = p(1:numStrt,:);
            % generates optimal PLHS based on 'maxmin' criterion.
            % by default the number of iteration is 10.
        else
            strtPoints = lhsdesign(numStrt, numDim,'criterion','maximin', 'iterations',100 );
        end
    case 'Sobolseq'
        p = sobolset(numDim,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
        p = scramble(p,'MatousekAffineOwen'); %MatousekAffineOwen scramling strategy. For no scrambling disable this line.
        strtPoints = p(1:numStrt,:);
    case 'Halton'
        p = haltonset(numDim,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
        p = scramble(p,'RR2'); %PR2 scramling strategy. For no scrambling disable this line.
        strtPoints = p(1:numStrt,:);
    otherwise
        strtPoints = lhsdesign(numStrt, numDim,'criterion','maximin', 'iterations',100 );
end
%% Calculate Elementary Effects (ElEf)
ElEf = zeros(numStrt, numDim);
fprintf('A Morris experiment started: number of restarts (chains) = %g, step size = %g, total number of model runs = %g. \n', numStrt, stpSize, MorrisCost);
for i = 1: numStrt
    fprintf('Chain #%g started. Running model %s %g times...', i, funFile, numDim + 1 );
    tic; [ ElEf(i, :) ] = MorrisCalc (strtPoints(i, :), stpSize, lb, ub, funPath, funFile); time = toc;
    fprintf(' Chain finished in %g seconds.\n', time);
end
Morris_out.ACE = mean(ElEf, 1);      % mean ACtual Elementary effect
Morris_out.SDE = std(ElEf, 1);       % standard deviation of elementary effect
Morris_out.ABE = mean(abs(ElEf), 1); % mean ABsolute Elementary effect
Morris_out.SQE = mean(ElEf .^ 2, 1); % mean SQuare Elementary effect
Morris_out.rnkACE = factorRanking(Morris_out.ACE);
Morris_out.rnkABE = factorRanking(Morris_out.ABE);
Morris_out.rnkSQE = factorRanking(Morris_out.SQE);
%% Bootstrap
if bootstrapFlag == 1;
    [ Morris_out.bootstrap ] = bootstrapMorris (ElEf, bootstrapSize, confdLvl, Morris_out.rnkACE, Morris_out.rnkABE, Morris_out.rnkSQE, numGrp);
end
%% Store Results
save ('Results_Morris', 'Morris_inp', 'Morris_out');
end
%% ************************************************************************
function [ EfEl ] = MorrisCalc (X0, stepSize, lb, ub, funPath, funFile)
currDir = pwd;
cd(funPath);
numDim = length(lb);
temp = ( [[1 : 1 : numDim ]'  rand(numDim, 1)] );
temp = sortrows(temp, 2);
shuffled_factors = temp(:, 1);
X_orgScale = X0 .* ( ub' - lb') + lb';
y0 = feval(funFile, X_orgScale );
for i = 1 : numDim
    signMultiplier = sign( rand() - 0.5 );
    if i == 1
        X(i, :) = X0;
        temp1 = X(i, shuffled_factors(i) ) + signMultiplier * stepSize;
        temp2 = X(i, shuffled_factors(i) ) - signMultiplier * stepSize;
        if temp1 <= 1 && temp1 >= 0
            X(i, shuffled_factors(i) ) = temp1;
        else
            X(i, shuffled_factors(i) ) = temp2;
        end
        X_orgScale = X(i, :) .* ( ub' - lb') + lb';
        y(i) = feval(funFile, X_orgScale );
        EfEl(shuffled_factors(i)) = ( y(i) - y0 );
    else
        X(i, :) = X(i - 1, :);
        temp1 = X(i, shuffled_factors(i) ) + signMultiplier * stepSize;
        temp2 = X(i, shuffled_factors(i) ) - signMultiplier * stepSize;
        if temp1 <= 1 && temp1 >= 0
            X(i, shuffled_factors(i) ) = temp1;
        else
            X(i, shuffled_factors(i) ) = temp2;
        end
        X_orgScale = X(i, :) .* ( ub' - lb') + lb';
        y(i) = feval(funFile, X_orgScale );
        EfEl(shuffled_factors(i)) = ( y(i) - y(i - 1) );
    end
end
cd (currDir);
end
%% ------------------------------------------------------------------------
% function [lb  ub] = readFactorSpace(funPath, facrSpcFile)
% currDir = pwd;
% cd(funPath);
% temp = importdata(facrSpcFile);
% cd(currDir);
% lb = temp.data(:, 2);
% ub = temp.data(:, 3);
% headers = temp.colheaders;
% end
%% ------------------------------------------------------------------------
function [ bootstrap ] = bootstrapMorris (ElEf, bootstrapSize, confdLvl, rnkBnchmrkACE, rnkBnchmrkABE, rnkBnchmrkSQE, numGrp)

[ N  numDim ] = size(ElEf);

[randstream1] = RandStream.create('mrg32k3a','NumStreams',1, 'Seed', 1234567);
for k = 1 : bootstrapSize
    rows = ceil ( rand(randstream1, N, 1) * N );
    ElEf_boot = ElEf(rows, :);
    
    ACE(k, :) = mean(ElEf_boot, 1);      % mean ACtual Elementary effect
    SDE(k, :) = std(ElEf_boot, 1);       % standard deviation of elementary effect
    ABE(k, :) = mean(abs(ElEf_boot), 1); % mean ABsolute Elementary effect
    SQE(k, :) = mean(ElEf_boot .^ 2, 1); % mean SQuare Elementary effect
    
    rnkACE(k, :) = factorRanking( ACE(k, :) );
    rnkSDE(k, :) = factorRanking( SDE(k, :) );
    rnkABE(k, :) = factorRanking( ABE(k, :) );
    rnkSQE(k, :) = factorRanking( SQE(k, :) );
end

ACE_sorted = sort(ACE, 1);
SDE_sorted = sort(SDE, 1);
ABE_sorted = sort(ABE, 1);
SQE_sorted = sort(SQE, 1);

bootstrap.ACE_low = ACE_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.ACE_upp = ACE_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );
bootstrap.SDE_low = SDE_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.SDE_upp = SDE_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );
bootstrap.ABE_low = ABE_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.ABE_upp = ABE_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );
bootstrap.SQE_low = SQE_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.SQE_upp = SQE_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );

for D = 1 : numDim
    bootstrap.rel_ACE(1, D) = length ( find(  rnkACE(:, D) == rnkBnchmrkACE(1, D) ) ) / bootstrapSize;
    bootstrap.rel_ABE(1, D) = length ( find(  rnkABE(:, D) == rnkBnchmrkABE(1, D) ) ) / bootstrapSize;
    bootstrap.rel_SQE(1, D) = length ( find(  rnkSQE(:, D) == rnkBnchmrkSQE(1, D) ) ) / bootstrapSize;
end
end
%% ************************************************************************
function [ rank ] = factorRanking(SAindices)
[sorted, order] = sort(SAindices, 'descend');
temp = [ order; 1: length(SAindices) ]';
temp2 = sortrows(temp, 1)';
rank = temp2(2, :);
end
%% ************************************************************************
function rank_grp = groupRanking ( rank_indvl, numGrp )
numDim = length (rank_indvl);
grpSize = round ( numDim / numGrp );
grpNum = 1; temp = 0;
for rankNum = 1 : numDim
    if temp == grpSize
        temp = 0;
        grpNum = grpNum + 1;
    end
    temp = temp + 1;
    rank_grp ( 1,  rank_indvl == rankNum  ) = grpNum;
end
end