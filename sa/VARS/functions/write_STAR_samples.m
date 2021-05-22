function write_STAR_samples(pivots, X_star, X_scaled_mat, VARS_inp )

if exist (VARS_inp.outFldr, 'dir') ~= 7
       mkdir ( VARS_inp.outFldr );
end
filename1 = strcat (VARS_inp.outFldr, '/STAR_out' );

save ( filename1, 'pivots', 'X_star', 'X_scaled_mat');

filename = strcat (VARS_inp.outFldr, '/STAR_out', '.smp' );
fid = fopen(filename,'w');

fprintf(fid,' -------- Star-based (STAR) Sampling Output File -------- \r\n');
fprintf(fid,' STAR settings: \r\n');

fprintf(fid,' Number of Star Points = %g \r\n', VARS_inp.numStars);
fprintf(fid,' minimum h = %g \r\n', VARS_inp.grdSize);
fprintf(fid,'\r\n        factor#1        factor#2        factor#3        -------> \r\n');

for i = 1: size(X_scaled_mat, 2)
    for j = 1 : size(X_scaled_mat{i}, 1)
        fprintf(fid,'%15g ', X_scaled_mat{i}(j, :));
        fprintf(fid,'\r\n');
    end
end

fclose(fid);

end