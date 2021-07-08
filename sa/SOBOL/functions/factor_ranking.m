function [ rank ] = factor_ranking(SAindices)
% Author: Saman Razavi in November 2019

[~, order] = sort(SAindices, 'descend');
temp = [ order; 1: length(SAindices) ]';
temp2 = sortrows(temp, 1)';
rank = temp2(2, :);
end