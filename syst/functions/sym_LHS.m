function [sample] = sym_LHS( n , d , criterion , maxEval)
%%
% =============================================================== %
% Description:
% This code generates Symmetric Latin Hypercube Design (SymLHS) of
% n datapoints in the d-dimensional hypercube [0,1].
% The algorithm is originally developed by Ye et al. (2000)
%----------------------------------------------------------------
% Programmed by Razi Sheikholeslami at GIWS, School of Environment
% and Sustainability, University of Saskatchewan

% Programming date: April 2018    version: 1.0
% E-mail: razi.sheikholeslami@usask.ca
% ----------------------------------------------------------------
% Original paper:
% 1) Ye, K.Q., Li, W., Sudjianto, A., 2000. Algorithmic construction of
%    optimal symmetric Latin hypercube designs.
%     Journal of Statistical Planning and Inference. 90 (1), 145–159.
% =============================================================== %
% Inputs:
%             n ==========>>>  number of sample points
%             d ==========>>>  number of parameters
%           maxEval ==========>>> maximum number of evaluations in
%           optimization process (default value is 50)
%          criterion ==========>>>  optimization criterion = 'maxmin' or
%          'correlation' (default is 'maxmin')
%  Outputs:
%             sample ==========>>>  n x d Symmetric Latin Hypercube Sample

% =============================================================== %
% Notes:
% 1)In this verversion in order to obtain the optimal SymLHS the maxmin
% and correlatoin criteria are considered i.e. the objective is maximizing
% the minimum point-to-point distance of a desing or minimizing the correlation.
% In optimization process a simple random search is utilized to find an
% optimal design.

%%  Check inputs

if ~isscalar(n); error('''n'' must be a scalar'); end
if n<=1; error('''n'' must be larger than 1' ); end
if abs(n-round(n)); error('''n'' must be integer'); end

if ~isscalar(d); error('''d'' must be a scalar'); end
if d<=1; error('''d'' must be larger than 1' ); end
if abs(d-round(d)); error('''d'' must be integer'); end

%% Default optional inputs

if nargin<4 && nargin>2
    maxEval = 50;
elseif  nargin<3
    maxEval = 50;
    criterion = 'maxmin';
end

%% Check the optimization criterion
switch (criterion)
    case 'maxmin'
        best_sample = generate_symLHS(n,d);
        best_sample_cost = get_min_dist(best_sample);    
        for itr = 2:maxEval
            new_sample = generate_symLHS(n,d);
            new_sample_cost = get_min_dist(new_sample);
            
            if new_sample_cost > best_sample_cost
                best_sample = new_sample;
                best_sample_cost = new_sample_cost;
            end
        end
        sample = best_sample;
        
    case 'correlation'
        best_sample = generate_symLHS(n,d);
        best_sample_cost = get_corr(best_sample);
        
        for itr = 2:maxEval
            new_sample = generate_symLHS(n,d);
            new_sample_cost = get_corr(new_sample);        
            if new_sample_cost < best_sample_cost
                best_sample = new_sample;
                best_sample_cost = new_sample_cost;
            end
        end
        sample = best_sample;
end
%% Main code
    function [sample] = generate_symLHS(n,d)
        X = zeros(n,d);
        if mod(n,2)==0
            start = 1;
        else
            start = 2;
            X(1,:) = (n+1)/2;
        end
        for i = start:2:n
            X(i,:) = rand_perm_interval(1,n,d);
            if i>1
                for j = 1:d
                    eq_indx = find(hist(X(1:i,j),unique(X(1:i,j)))>1);
                    while eq_indx >= 1
                        X(i,j) = rand_perm_interval(1,n,1);
                        eq_indx = find(hist(X(1:i,j),unique(X(1:i,j)))>1);
                    end
                end
            end
            X(i+1,:) = n+1-X(i,:);
        end      
        X = unifrnd(X-1,X);
        sample = X/n;
    end
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
end