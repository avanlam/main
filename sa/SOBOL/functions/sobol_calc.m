function [ FO, TO, Vy, MUy ] = sobol_calc (yA, yB, yC)
% Author: Saman Razavi in November 2019

[ N,  numDim ] = size(yC);
f0 = mean([ yA; yB ]);  % sample mean - approximation of the mean of response surface
MUy = f0;
Vy = ( [ yA; yB ]' * [ yA; yB ] ) / ( 2 * N ) - f0 ^ 2; % sample variance - approximation of the variance of response surface
for i = 1 : numDim
    FO(i) = ( ( (yA-f0)' * (yC(:, i) -f0)) / N ) / Vy;     % First-Order effect
    TO(i) = 1 - ( ( yB' * yC(:, i) ) / N - f0 ^ 2 ) / Vy;  % Total-Order effect
end
end