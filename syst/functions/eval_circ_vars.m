function [out] = eval_circ_vars(factors)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

num_n = 100;
output = 3;
factors = [num_n, factors, output];

out = eval_circ(factors);

end

