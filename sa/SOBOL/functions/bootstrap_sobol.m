function [ bootstrap ] = bootstrap_sobol (yA, yB, yC, bootstrapSize, confdLvl, rnkFOBnchmrk, rnkTOBnchmrk, numGrp)

[ N,  numDim ] = size(yC);
[randstream1] = RandStream.create('mrg32k3a','NumStreams',1, 'Seed', 1234567);
for k = 1 : bootstrapSize
    rows = ceil ( rand(randstream1, N, 1) * N );
    yAboot = yA(rows);
    yBboot = yB(rows);
    yCboot = yC(rows, :);
    [ FO(k, :), TO(k, :), V(k, :), MU(k, :) ] = sobol_calc (yAboot, yBboot, yCboot);
     rnkFO(k, :) = factor_ranking(FO(k, :));
     rnkTO(k, :) = factor_ranking(TO(k, :));
end
FO_sorted = sort(FO, 1);
TO_sorted = sort(TO, 1);

bootstrap.FO_low = FO_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.FO_upp = FO_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );

bootstrap.TO_low = TO_sorted ( round ( bootstrapSize * ( ( 1 - confdLvl ) / 2 ) ) , : );
bootstrap.TO_upp = TO_sorted ( round ( bootstrapSize * ( 1 - ( ( 1 - confdLvl ) / 2 ) ) ) , : );

for D = 1 : numDim
    bootstrap.rel_FO(1, D) = length ( find( rnkFO(:, D) == rnkFOBnchmrk(1, D) ) ) / bootstrapSize;
    bootstrap.rel_TO(1, D) = length ( find( rnkTO(:, D) == rnkTOBnchmrk(1, D) ) ) / bootstrapSize;
end
% ************ based on group ranking 
rank_benchmark_grp = group_ranking ( rnkTOBnchmrk , numGrp );
for iter = 1 : bootstrapSize
    rnkTO_grp( iter, :) = group_ranking ( rnkTO(iter, :) , numGrp );
end
for D = 1 : numDim
    relGrp(1, D) = length ( find( rnkTO_grp(:, D) == rank_benchmark_grp(1, D) ) ) / bootstrapSize;
end
 
end