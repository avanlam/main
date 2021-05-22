function [ yA, yB, yC ] = func_eval (A, B, C, funPath, funFile, numN, chosen_o)
currDir = pwd;
cd(funPath);
[ N , numDim ] = size(A);
SobolCost = N * (numDim + 2); % total number of function evaluations (model runs)

yA = zeros(N, 1);
yB = zeros(N, 1);
yC = zeros(N, numDim);
fprintf(['A Sobol experiment started for %g neuron(s): size of base sample = %g, total number of model runs = %g. \n'], numN, N, SobolCost );
for j = 1 : N
    fprintf('Group run (base sample) #%g started. Running model %s %g times...', j, funFile, numDim + 2 );
    tic;
    yA(j, 1) = feval(funFile, [numN, A(j, :), chosen_o] );
    yB(j, 1) = feval(funFile, [numN, B(j, :), chosen_o] );
    for i = 1 : numDim
        yC(j, i) = feval(funFile, [numN, C{i}(j, :), chosen_o] );
    end
    time = toc;
    fprintf(' Group run finished in %g seconds.\n', time);
end
cd (currDir);
end