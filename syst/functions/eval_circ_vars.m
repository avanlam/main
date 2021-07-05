function [out] = eval_circ_vars(factors)
% EVAL_CIRC_VARS : Adaptation before EVAL_CIRC to incorporate the number of
% oscillators and the chosen synchronisation output

num_n = 100;
output = 3;     % Kuramoto order parameter
factors = [num_n, factors, output];

out = eval_circ(factors);

end

