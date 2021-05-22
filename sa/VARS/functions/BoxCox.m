function [trans_data, lambda_opt] = BoxCox( x )

% =============================================================== %
% Description:
% Box-Cox transformation
%
% By selection of an appropriate value for lambda, the data x is transformed to T
% which has a distribution that is closerto being symmetrical (i.e., having
% skewness close to zero). This is easily done by finding the value of lambda that
% maximizes a log-likelihood function using an optimization approach (Sakia, 1992;
% Box and Cox, 1964). Here, we used the Box-Cox transformation, enabled
% by a Nelder-Mead simplexdirect search optimization algorithm of MATLAB, to normalize
% the distribution of the sensitivity indices.

%----------------------------------------------------------------
% Programmed by Razi Sheikholeslami at GIWS, School of Environment
% and Sustainability, University of Saskatchewan, 2018.
%
% E-mail: razi.sheikholeslami@usask.ca
% ----------------------------------------------------------------
% Ref:
% 1)Box, G.E. and Cox, D.R., 1964. An analysis of transformations. 
%   Journal of the Royal Statistical Society. Series B (Methodological), 26(2), 211–252.
% 2)Sakia, R.M., 1992. The Box-Cox transformation technique: a review. 
%   The Statistician, 41(2), 169–178. https://doi.org/10.2307/2348250
% =============================================================== %

% Inputs
%            x ==========>>> data vector (n X 1)

%  Outputs:
%             trans_data ==========>>>  normaized data vector (n X 1)
%             lambda_opt ==========>>>  optimum lambda
% =============================================================== %

%%

% "fminsearch" algorithm is used here:
%Choice of the accuracies
precf = 10^-3;%'Objective' accuracy
precx = 10^-3;%'Simplex' accuracy
options = optimset('TolF',precf,'TolX',precx,'MaxFunEvals',...
    2000,'MaxIter',2000,'Display','off');
[lambda_opt] = sub_fun1( x, options );

% Transforming data using an optimum lambda
z1 = find(lambda_opt ~= 0);
z2 = find(lambda_opt == 0);
new_x = x * ones(1, length(lambda_opt));
new_lambda2 = (lambda_opt * ones(length(x), 1)')';

trans_data(:, z1) = ((new_x(:, z1).^new_lambda2(:, z1))-1) ./ ...
    new_lambda2(:, z1);
trans_data(:, z2) = log(new_x(:, z2));

end

function [best_lambda] = sub_fun1( x, options )

best_lambda = fminsearch(@loglike, 0, options);

    function Fobj = loglike(lambda)
        
        % Transforming data using a given lambda
        nz1 = find(lambda ~= 0);
        nz2 = find(lambda == 0);
        new_x = x * ones(1, length(lambda));
        new_lambda = (lambda * ones(length(x), 1)')';
        
        x_transformed(:, nz1) = ((new_x(:, nz1).^new_lambda(:, nz1))-1) ./ ...
            new_lambda(:, nz1);
        x_transformed(:, nz2) = log(new_x(:, nz2));
        
        %Find the lambda that minimizes the Log-Likelihood function
        Fobj = -(-((length(x))/2).*log(std(x_transformed', 1, 2).^2) + (lambda-1)*(sum(log(x))));

    end
end
