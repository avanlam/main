function [opt_PLHS,slices] = PLHSdesign(n,d,t,maxIter,criterion)
%% =============================================================== %
% Description:
% Generate an Optimal Progressive Latin Hypercube Sample (PLHS)
%
% The proposed PLHS is composed of a series of smaller slices generated in
% a way that the union of these slices from the beginning to the current stage
% optimally preserves the intended distributional properties and at the same
% time achieves maximum space-filling.

% The algorithm is proposed by Sheikholeslami & Razavi (2017)
%----------------------------------------------------------------
% Programmed by Razi Sheikholeslami at GIWS, School of Environment
% and Sustainability, University of Saskatchewan
%
% E-mail: razi.sheikholeslami@usask.ca
% ----------------------------------------------------------------
% Original paper:
% 1) R. Sheikholeslami & S. Razavi (2017):
% Progressive Latin Hypercube Sampling: An Efficient Approach for Robust
% Sampling-based Analysis of Environmental Models,Env. Model. Software.
% =============================================================== %
% Notes:
% 1)In this verversion in order to find an optimal PLHS
%   a random search optimization method is utilized.
% =============================================================== %
% Inputs:
%            n ==========>>> number of sample points
%            d ==========>>> number of parameters
%            t ==========>>> number of slices/sub-samples (n = m*t)
%      maxIter ==========>>> maximum number of iterations in random search
%    criterion ==========>>> optimization criterion = 'maxmin' or 'correlation'

%  Outputs:
%             opt_PLHS ==========>>>  Optimal PLHS
%             slices ==========>>>  Slices/sub-samples
% =============================================================== %

%%%%%%%%%%%%%%
%% Check inputs
%%%%%%%%%%%%%%

if ~isscalar(n); error('''n'' must be a scalar'); end
if n<=1; error('''n'' must be larger than 1' ); end
if abs(n-round(n)); error('''n'' must be integer'); end

if ~isscalar(d); error('''d'' must be a scalar'); end
if d<=1; error('''d'' must be larger than 1' ); end
if abs(d-round(d)); error('''d'' must be integer'); end
if ~mode(n,t);  error('''the remainder after division of n by t must be zero'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5
    criterion = 'maxmin';
end
if nargin<4
    maxIter = 50;
end
if nargin<3
    error('not enough input params');
end

wrrt=0; %optional. 1 to write the samples in text files; 0 do not
%% <><><><><><><><><><><><><><><><><>main part<><><><><><><><><><><><><><<><>

m=n/t;
for ii=1:maxIter
    %generating optimal SLHS
    [X{ii}] = SLHSdesign(n,d,t,2*maxIter,criterion);
    %generating PLHS using a greedy search
    [opt_s{ii},~,F(1,ii),F(2,ii)] = Greedy_Plhs(n,t,X{ii});
end
%finding the optimal PLHS
cost = F(2,:);
[~,index] = min(cost);

opt_PLHS = opt_s{index};
slices = mat2cell(opt_PLHS,m*ones(1,n/m));

%writing the samples in  text files
if wrrt==1
    fid = fopen('SampleOut.txt','wt');
    for kk = 1:size(opt_PLHS,1)
        fprintf(fid,'%f\t',opt_PLHS(kk,:));
        fprintf(fid,'\n');
    end
    fclose(fid);
    
    for wrt = 1:t
        filename = strcat ('slice#', num2str(wrt),'.txt');
        fid = fopen(filename,'wt');
        for jj = 1:m
            fprintf(fid,'%f\t',slices{wrt,1}(jj,:));
            fprintf(fid,'\n');
        end
        fclose(fid);
    end
end

end

%%
function [PLHS, PLHS_slice,F1, F2] = Greedy_Plhs (n,t,Sample)

% =============================================================== %
% Description:
% Generate a Progressive Latin Hypercube Sample (PLHS) from an Optimal
% Sliced Latin Hypercube Design (SLHSdesign)using a greedy algorithm.
%
% The algorithm is proposed by Sheikholeslami & Razavi (2017)
%----------------------------------------------------------------
% Programmed by Razi Sheikholeslami at GIWS, School of Environment
% and Sustainability, University of Saskatchewan
%
% E-mail: razi.sheikholeslami@usask.ca
% ----------------------------------------------------------------
% Original paper:
% 1) R. Sheikholeslami & S. Razavi (2017):
% Progressive Latin Hypercube Sampling: An Efficient Approach for Robust
% Sampling-based Analysis of Environmental Models,Env. Model. Software.
% =============================================================== %
% Notes:
% 1)In this verversion in order to find an optimal ordering of slices
%   the Nearest neighbor search (NNS)method is utilized.
% =============================================================== %
% Inputs:
%             n ==========>>>  number of sample points
%             t ==========>>>  number of slices/sub-samples (n = m*t)
%           Sample ==========>>> n*d SLHSdesign
%  Outputs:
%             PLHS ==========>>>  Progressive Latin Hypercube Sample
%             PLHS_slice ==========>>>  Slices/sub-samples

%             F1 =============>>> objcetive function value befor
%             optimization
%             F2 =============>>> objcetive function value after
%             optimization
%% =============================================================== %

m = n/t;
sample = Sample;
sub_sample = mat2cell(sample ,m*ones(1,t));

stage_p = [];
for i = 1:t
    stage_p = [stage_p ; sub_sample{i,1}];
    f(i) = LHD_cost(stage_p);
end
F1 = mean(f);

k = 1;
edge = ones(t,t);
while k<t
    for j = k+1:t
        stage = [sub_sample{k,1} ; sub_sample{j,1}];
        edge(k,j) = LHD_cost(stage);
    end
    k = k+1;
end
[~, minIndex] = min(edge(:));
[q(1), q(2)] = ind2sub(size(edge), minIndex);
slice = [sub_sample{q(1),1} ; sub_sample{q(2),1}];
X = sub_sample;
X{q(1),1}=[];
X{q(2),1}=[];
Y = cell2mat(X);
X = mat2cell(Y,m*ones(1,t-2));

t = t-2;
count = 3;
while t>0
    dist=[];
    for j = 1:t
        stage = [slice ; X{j,1}];
        dist(j) = LHD_cost(stage);
    end
    [~,q(count)]=min(dist);
    slice = [slice ; X{q(count),1}];
    X{q(count),1}=[];
    Y = cell2mat(X);
    X = mat2cell(Y,m*ones(1,t-1));
    t = t-1;
    count = count+1;
    
end

PLHS = slice;
PLHS_slice = mat2cell(PLHS,m*ones(1,n/m));

p_stage = [];
for i = 1:(n/m)
    p_stage = [p_stage ; PLHS_slice{i,1}];
    ff(i) = LHD_cost(p_stage);
end

F2 = mean(ff);
end

%% <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
function F = LHD_cost(X) %X: N rows, d col
[N, d]= size(X);
v = 1:N;
edges = 0:1/N:1;
Y = discretize(X',edges);
f=zeros(1,d);

for i = 1:d
    u = unique(Y(i,:));
    p = hist(u,v);
    f(i) = sum(p);
end
F = -sum(f);
end
%%
function [sample,sub_samples] = SLHSdesign(n,d,t,maxEval,criteria)

% =============================================================== %
% Description:
% Generate an Optimal Sliced Latin Hypercube Design (SLHSdesign) of
% n datapoints in the d-dimensional hypercube [0,1] which is a union of
% t small Latin hypercubes with m=n/t sample points.
% is an efficient
% The algorithm is originally developed by Shan Ba et al. (2014)
%----------------------------------------------------------------
% Programmed by Razi Sheikholeslami at GIWS, School of Environment
% and Sustainability, University of Saskatchewan

% Programming date: Sep 2015    version: 1.0
%                   Oct 2015    version: 1.1
% E-mail: razi.sheikholeslami@usask.ca
% ----------------------------------------------------------------
% Original paper:
% 1) Shan Ba, William A. Brenneman & William R. Myers (2014):
% Optimal Sliced Latin Hypercube Designs, Technometrics,
% DOI: 10.1080/00401706.2014.957867
% =============================================================== %
% Notes:
% 1)In this verversion in order to obtain the optimal design the maxmin
% and correlatoin criteria are considered i.e. the objective is maximizing
% the minimum point-to-point distance of a desing or minimizing the correlation.
% In optimization process a simple random search is utilized to find an
% optimal design.
%
% =============================================================== %
% Inputs:
%             n ==========>>>  number of sample points
%             d ==========>>>  number of parameters
%             t ==========>>>  number of slices/sub-samples (n = m*t)
%           maxEval ==========>>> maximum number of evaluations in optimization process
%          criteria ==========>>>  optimization criteria = 'maxmin' or 'correlation'

%  Outputs:
%             sample ==========>>>  Optimal Sliced Latin Hypercube Sample
%             sub_sample ==========>>>  Slices/sub-samples
% =============================================================== %

%%%%%%%%%%%%%%
%% Check inputs
%%%%%%%%%%%%%%

if ~isscalar(n); error('''n'' must be a scalar'); end
if n<=1; error('''n'' must be larger than 1' ); end
if abs(n-round(n)); error('''n'' must be integer'); end

if ~isscalar(d); error('''d'' must be a scalar'); end
if d<=1; error('''d'' must be larger than 1' ); end
if abs(d-round(d)); error('''d'' must be integer'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Default optional inputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin<5
    maxEval = 50;
    criteria = 'maxmin';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Check the optimization criterion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = n/t;
switch (criteria)
    case 'maxmin'
        best_sample = generate_sample(n,d,t,m);
        best_sub_sample = mat2cell(best_sample ,m*ones(1,t));
        best_sample_cost = get_min_dist(best_sample);
        best_sub_sample_cost = get_subdist(best_sub_sample);
        Cost = (best_sample_cost + best_sub_sample_cost)/2;
        
        for itr = 2:maxEval
            new_sample = generate_sample(n,d,t,m);
            new_sub_sample = mat2cell(new_sample,m*ones(1,t));
            new_sample_cost = get_min_dist(new_sample);
            new_sub_sample_cost = get_subdist(new_sub_sample);
            new_Cost = (new_sample_cost + new_sub_sample_cost)/2;
            
            if new_Cost > Cost
                best_sample = new_sample;
                Cost = new_Cost;
            end
        end
        sample = best_sample;
        
    case 'correlation'
        best_sample = generate_sample(n,d,t,m);
        best_sub_sample = mat2cell(best_sample ,m*ones(1,t));
        best_sample_cost = get_corr(best_sample);
        best_sub_sample_cost = get_subcorr(best_sub_sample);
        Cost = (best_sample_cost + best_sub_sample_cost)/2;
        
        for itr = 2:maxEval
            new_sample = generate_sample(n,d,t,m);
            new_sub_sample = mat2cell(new_sample ,m*ones(1,t));
            new_sample_cost = get_corr(new_sample);
            new_sub_sample_cost = get_subcorr(new_sub_sample);
            new_Cost = (new_sample_cost + new_sub_sample_cost)/2;
            if new_Cost < Cost
                best_sample = new_sample;
                Cost = new_Cost;
            end
        end
        sample = best_sample;
end

Y=sample;
sub_samples = mat2cell(Y,m*ones(1,t));
%sub_sample=cell(t,1);
%for j = 1:t
%   sub_sample{j,1} = Y{j,1}';
%end
end

%<><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
%%
function perm = rand_perm_interval(lb,ub,t)
size = ub-lb+1;
perm = lb:ub;
for k = 2:size
    index1 = ceil(rand*k);
    index2 = perm(k);
    perm(k) = perm(index1);
    perm(index1) = index2;
end
perm = perm(1:t);
end
%%
function X = generate_sample(n,d,t,m)

Y = cell(1,t);
q = zeros(m,n);
q = logical(q);

for i = 1:t
    for j = 1:d
        Y{i}(j,:) = randperm(m);
    end
end

X = cell2mat(Y);

for j = 1:d
    for k = 1:m
        q(k,:) = (X(j,:)==k);
    end
    for kk = 1:m
        lb = (kk-1)*t+1;
        ub = kk*t;
        p = rand_perm_interval(lb,ub,t);
        X(j,q(kk,:)) = p;
    end
end
X = unifrnd(X-1,X);
X = X/n;
X = X';
end
%%
function mean_min_dist =get_subdist(X)
min_dist =zeros(1,size(X,1));
for i=1:size(X,1)
    min_dist(i) = get_min_dist(X{i,1});
end
mean_min_dist = mean(min_dist);
end
%%
function mean_corr =get_subcorr(X)
min_corr =zeros(1,size(X,1));
for i=1:size(X,1)
    min_corr(i) = get_corr(X{i,1});
end
mean_corr = mean(min_corr);
end
%%
function min_dist = get_min_dist(X)
% Maximize the minimum point-to-point difference

[~,dist] = knnsearch(X,X,'k',3);
min_dist = min(dist(:,2));
%min_dist = min(pdist(X,'correlation'));

end
%%
function min_corr = get_corr(X)
% Minimize the sum of between-column squared correlations

R = corrcoef(X);
min_corr = sum(sum(triu(R,1).^2));
end
