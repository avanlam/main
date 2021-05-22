function [RSA_inp, RSA_out] = main_RSA( )
%% =============================================================== %
% Description:
% This function is a part of the VARS-TOOL software. It performs Sensitivity Analysis based on
% Regional Sensitivity Analysis algorithm enabled by bootsrap strategy and plots the corresponding
% behavioral and non-behavioral sample sets.
% =============================================================== %
%  Outputs:
%         RSA_out.SA_index: sensitivity indices calculated using the Kolmogorov-Smirnov (K-S) measure
%         RSA_out.rnk:  factor importance ranking (most important factor is the 1st rank)
%         RSA_out.bootstrap: results of bootstrapping
% =============================================================== %
% Note:
% The K-S measure calculates the maximum vertical distance between the two CDFs.
% =============================================================== %
% References:
% (1) Hornberger, G., and Spear, R. (1981). “Approach to the preliminary
%     analysis of environmental systems.” J. Environ. Manage., 12, 7–18.
% (2) Stephens, M. A. (1992). “Introduction to Kolmogorov (1933) on the empirical determination
%     of a distribution.” Breakthroughs in statistics. Springer series in statistics
%     (perspectives in statistics), Kotz, G. and Johnson, N. L., eds., Springer, New York.
%----------------------------------------------------------------
% Programmed by: Razi Sheikholeslami, April 2018.
% Global Institute for Water Security,
% School of Environment and Sustainability,
% University of Saskatchewan
% E-mail: razi.sheikholeslami@usask.ca
% ----------------------------------------------------------------
%% Control Parameters and Specifications
N = 1000;                             % 1. Sample size: total number of function evaluation will be N
seedNum = 123456789;                  % 2. Seed for random number generator; put [] for autmatic randomization
funPath = 'C:\2019_11_16_VARS-Tool-v2.1\HBV-SASK'; % 3. Folder address: that includes model file, factor space
funFile = 'eval_HBV_SASK';            % 4. Model/function file: MATLAB m-file without .m extension
smplMtd = 'LHS';                      % 5. Sampling Method: RND, LHS, SymLHS, PLHS, SobolSeq, or Halton for generation of star centers;                                       %    if blank, default is LHS
sliceSize = 10;                       % 6. Slice Size: needed only for PLHS sampling strategy
bootstrapFlag = 1;                    % 7. Bootstrap Flag: enter 1 to bootstrap, 0 not to bootstrap
DoPlot = 1;                           % 8. plotting flag, 1 = do plot.
bootstrapSize = 1000;                 % 9. Bootstrap Size: number of sampling iterations with replacement              
                                      %    if bootstrap flag = 0, this will be ignored.
confdLvl = 0.9;                       % 10. Confidence Level: to report bootstrap-based lower and upper bounds
                                      %   on VARS results; if bootstrapFlag = 0, this line will be ignored.
%% Store RSA Control Parameters and Specifications
RSA_inp.N = N;
RSA_inp.seedNum = seedNum;
RSA_inp.funFile = funFile;
RSA_inp.funPath = funPath;
RSA_inp.smplMtd = smplMtd;
%% Randomization
if isempty( seedNum ) == false
    rand('state',seedNum);
end
%% Read variables
factors = read_factorSpace(funPath); % read in the factor space
lb = factors.lb; ub = factors.ub;
num_param = length(lb);  % number of factors

%% Generate the base sample from a unit hypercube
switch smplMtd
    case 'RND'
        sample = rand(N, num_param); %uniformly distributed random numbers
    case 'LHS'
        sample = lhsdesign(N, num_param,'criterion','maximin', 'iterations',100 );
    case 'SymLHS'
        sample = SymLHS( N , num_param,'maxmin', 10 );
    case 'PLHS'
        numSlices = ceil(N/sliceSize);
        if numSlices > 1
            smpSize = numSlices * sliceSize;
            
            p = PLHSdesign(smpSize, num_param, numSlices, 10, 'maxmin');
            sample = p(1:N,:);
            % generates optimal PLHS based on 'maxmin' criterion.
            % by default the number of iteration is 10.
        else
            sample = lhsdesign(N, num_param,'criterion','maximin', 'iterations',50 );
        end
    case 'Sobolseq'
        p = sobolset(num_param,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
        p = scramble(p,'MatousekAffineOwen'); %MatousekAffineOwen scramling strategy. For no scrambling disable this line.
        sample = p(1:N,:);
    case 'Halton'
        p = haltonset(num_param,'Skip',1e3,'Leap',1e2); %skip the first 1000 values, and then retain every 101st point
        p = scramble(p,'RR2'); %PR2 scramling strategy. For no scrambling disable this line.
        sample = p(1:N,:);
    otherwise
        sample = lhsdesign(N, num_param,'criterion','maximin', 'iterations',100 );
end

%% Scale sample set onto the actual factor space
param_set = sample(:, 1 : num_param);
for j = 1 : N
    param_set(j, :) = param_set(j, :) .* ( ub' - lb') + lb';
end
%% Run the function/model for all points in param_set and return the function/model response
[model_output] = funcEval(param_set, funPath, funFile);
[~,P]= size(model_output) ;

threshold = prctile(model_output, 95);  % A predefined threshold value for separating the model outputs into behavioral and nonbehavioral sets
% The default value is the 95th percentile of the model outputs
%% Check inputs
if size(param_set,1)-size(model_output,1)~=0; error('number of model runs should be equal to the sample size'); end
%% Calculate RSA indices and ranking
[SA_index, ranks] = calcRSA_main(model_output, param_set, N, threshold, DoPlot, P, factors);
RSA_out.SA_index = SA_index;
RSA_out.rnk = ranks;
%% Bootstrap
if bootstrapFlag == 1;
    [ RSA_out.bootstrap ] = bootstrapRSA (model_output, param_set, bootstrapSize, confdLvl, RSA_out.rnk,  threshold, P);
end
%% Store Results
save ('Results_RSA', 'RSA_inp', 'RSA_out');

end
%% Main code
function [SA_index, ranks] = calcRSA_main(model_output, param_set, N, threshold, DoPlot, P, factors)
% Finding behavioral and non-behavioral sets
num_param = size(param_set,2);
b_index = sum(model_output>repmat(threshold,N,1),2)==P;
param_b  = param_set(b_index,:);
b_length  = size(param_b,1);
param_nb = param_set(~b_index,:);
nb_length = size(param_nb,1);

SA_index = nan(1,num_param);
if b_length<=0 || nb_length<=0
    warning('Change the threshold value. Cannot assign any output given the current threshold value.')
else
    for dim = 1:num_param
        [Fbcdf, xbcdf]  = ecdf(param_b(:,dim)); % empirical CDF
        [Fnbcdf, xnbcdf] = ecdf(param_nb(:,dim)); % empirical CDF
        SA_index(dim)= kolmogorov_smirnov_distance(param_b(:,dim), param_nb(:,dim)); % Measure K-S distance between CDFs
        if DoPlot==1
            plotCDFs(num_param, xbcdf,Fbcdf, xnbcdf, Fnbcdf, dim, factors)
        end
    end
end
ranks = factorRanking(SA_index); %ranking the input factors

end

%% ************************************************************************
function [ model_output ] = funcEval (param_set, funPath, funFile)
currDir = pwd;
cd(funPath);
[ N , ~ ] = size(param_set);
RSACost = N; % total number of function evaluations (model runs)

model_output = zeros(N, 1);
fprintf('An RSA experiment started: size of base sample = %g, total number of model runs = %g. \n', N, RSACost );
for j = 1 : N
    fprintf('Model run #%g started. Running model %s %g ', j, funFile);
    tic;
    model_output(j, 1) = feval(funFile, param_set(j,:) );
    time = toc;
    fprintf(' Model run finished in %g seconds.\n', time);
end
cd (currDir);
end

%% Factor ranking
function [ rank ] = factorRanking(SAindices)
[~, order] = sort(SAindices, 'descend');
temp = [ order; 1: length(SAindices) ]';
temp2 = sortrows(temp, 1)';
rank = temp2(2, :);

end

%% Calculating K-S distance
function [KS, dist, nx] = kolmogorov_smirnov_distance(CDF1,CDF2)

l1 = length(CDF1);
l2 = length(CDF2);
dist = UniqueValues([CDF1;CDF2],1);
nx = length(dist);
dist = [dist, zeros(nx,2)];
incr = 1./l1;
for ii = 1:l1
    kk = find(dist(:,1)==CDF1(ii));
    dist(kk,2) = dist(kk,2)+incr;
end;

incr = 1./l2;
for ii = 1:l2
    kk = find(dist(:,1)==CDF2(ii));
    dist(kk,3) = dist(kk,3)+incr;
end;
dist(:,2) = cumsum(dist(:,2));
dist(:,3) = cumsum(dist(:,3));
diff = dist(:,2)-dist(:,3);

KS = max(abs(diff));          %K-S
end

function [uV,Freq,inds] = UniqueValues(VECTOR,FLAG)
if (nargin < 2), FLAG = []; end

getIndices = false;
if (nargout > 2), getIndices = true; end
if (isempty(FLAG)), FLAG = false; end

if (isempty(VECTOR))
    uV = [];
    Freq = [];
    inds = [];
    return
end

tol = eps*10.^6;
VECTOR = VECTOR(:);

if (getIndices)
    ind = (1:length(VECTOR))';
end

if (any(~isfinite(VECTOR)))  % remove NaN and infinite values
    i = find(~isfinite(VECTOR));
    VECTOR(i) = [];
    if (getIndices)
        ind(i) = [];
    end
end

uV = [];
Freq = [];

for s = 1:length(VECTOR)
    b = (abs(uV-VECTOR(s)) < tol);
    if (sum(b) > 0)
        Freq(b) = Freq(b) + 1;
    else
        uV = [uV; VECTOR(s)];
        Freq =  [Freq; 1];
    end
end

if (FLAG)
    [uV,s] = sort(uV);
    Freq = Freq(s);
end

if (getIndices)
    nval = length(uV);
    inds = zeros(nval,1);
    for v = 1:nval
        s = find(VECTOR == uV(v));
        inds(v) = ind(s(1));
    end
end
end
%% Bootstrap
function [ bootstrap ] = bootstrapRSA (model_output, param_set, bootstrapSize, confdLvl, rnks, threshold, P)

[N,~] = size(model_output);
[randstream1] = RandStream.create('mrg32k3a','NumStreams',1, 'Seed', 1234567);
SA_indexall = nan(bootstrapSize, size(param_set,2));

if P == 1 % print only in case of single-output models
    fprintf('Bootstrapping: Resampling Time (of %g Times) = ', bootstrapSize);
end

for k = 1 : bootstrapSize
    
    %%% to print on screen
    if P == 1 % print only in case of single-output models
        if k == 1
            fprintf('%g', k);
        elseif ( round(k/10) == k/10) || k == bootstrapSize
            for i = 0 : log10 ( k - 1 ); fprintf('\b'); end % delete previous counter display
            fprintf('%g', k);
        end
        if k == bootstrapSize; fprintf('\n'); end;
    end
    %%% ------------------
    
    rows = ceil ( rand(randstream1, N, 1) * N );
    yboot = model_output(rows, :);
    param_setboot = param_set(rows, :);
    [SA_indexall(k, :), rankall(k, :)] = calcRSA_main(yboot, param_setboot, N, threshold, 0, P, []);
end

idx = ~isnan(SA_indexall(:,1)) ;
if sum(idx)<bootstrapSize
    warning('Bootstrap-based SA indices were calculated using %d bootstrap replicates not %d replicates',sum(idx),bootstrapSize)
end
% Lower & Upper bound
SA_index_sorted = sort(SA_indexall(idx,:),1) ;
bootstrap.SA_index_low = SA_index_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.SA_index_upp = SA_index_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );

% Reliability
for D = 1 : size(param_set,2)
    bootstrap.rel_RSA(1, D) = length ( find( rankall(:, D) == rnks(1, D) ) ) / bootstrapSize;
end
end

%% Plotting function
function plotCDFs(num_param, xbcdf, Fbcdf, xnbcdf, Fnbcdf, dim, factors)

figure_column = 5; figure_column = min(floor(figure_column),num_param)  ;
figure_row = ceil(num_param/figure_column);

subplot(figure_row,figure_column,dim)
hold on
%ecdf(param_b(:,dim)); hold on; ecdf(param_nb(:,dim));
plot(xbcdf,Fbcdf,'LineWidth',1.5, 'LineStyle','-.'); hold on; plot(xnbcdf,Fnbcdf,'LineWidth',1.5)
ax=gca;
xlabel(factors.name{dim},'FontName','Times','FontAngle','italic','FontWeight','bold','FontSize',12)
ylabel('CDF','FontName','Times','FontWeight','bold','FontSize',12)
set(ax,'LineWidth',1.25,'FontSize',12);

title ={'Behavioural','Non-behavioural'};
if dim==1;
    if ~isempty(title); legend(title); end
end
box on
end
