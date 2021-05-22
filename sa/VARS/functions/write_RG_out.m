function write_RG_out(VARS_inp, starNum, RG)

if exist (VARS_inp.outFldr, 'dir') ~= 7
    mkdir ( VARS_inp.outFldr );
end

filename = strcat (VARS_inp.outFldr, '\VARS_RG_out_', num2str(starNum), '.txt' );

fid = fopen(filename,'w');

fprintf(fid,' --- Variogram Analysis of Response Surfaces (VARS) for Dynamical Systems Models Output File --- \r\n');
fprintf(fid,' VARS settings: \r\n');

fprintf(fid,' Target Number of Stars = %g \r\n', VARS_inp.numStars);
fprintf(fid,' Current Number of Stars = %g \r\n', starNum);

fprintf(fid,' minimum h = %g \r\n', VARS_inp.grdSize);

% print IVARS
fprintf(fid,'\r\n ***  IVARS Time Aggregate: \r\n');
fprintf(fid,'         H        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, VARS_inp.IVARS, RG.IVARS_tAgg);

fprintf(fid,'\r\n ***  Factor Rankings based on IVARS Time Aggregate: \r\n');
fprintf(fid,'         H        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, VARS_inp.IVARS, RG.IVARS_tAgg_rank);

fprintf(fid,'\r\n ***  IVARS Normalized Time Aggregate: \r\n');
fprintf(fid,'         H        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, VARS_inp.IVARS, RG.IVARS_ntAgg);

fprintf(fid,'\r\n ***  Factor Rankings based on IVARS Normalized Time Aggregate: \r\n');
fprintf(fid,'         H        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, VARS_inp.IVARS, RG.IVARS_ntAgg_rank);

% print VARS-TO
fprintf(fid,'\r\n ***  VARS-TO Time Aggregate: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
write2(fid, RG.TO_tAgg);

fprintf(fid,'\r\n ***  Factor Rankings based on VARS-TO Time Aggregate: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
write2(fid, RG.TO_tAgg_rank);

fprintf(fid,'\r\n ***  VARS-TO Normalized Time Aggregate: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
write2(fid, RG.TO_ntAgg);

fprintf(fid,'\r\n ***  Factor Rankings based on VARS-TO Normalized Time Aggregate: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
write2(fid, RG.TO_ntAgg_rank);

% print VARS-ABE
fprintf(fid,'\r\n ***  VARS-ABE Time Aggregate: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
write2(fid, RG.ABE_tAgg);

fprintf(fid,'\r\n ***  Factor Rankings based on VARS-ABE Time Aggregate: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
write2(fid, RG.ABE_tAgg_rank);

fprintf(fid,'\r\n ***  VARS-ABE Normalized Time Aggregate: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
write2(fid, RG.ABE_ntAgg);

fprintf(fid,'\r\n ***  Factor Rankings based on VARS-ABE Normalized Time Aggregate: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
write2(fid, RG.ABE_ntAgg_rank);

fclose(fid);
end
%% -----------------------------------------------
function write1(fid, rowTitle, matrix)
    for i = 1 : size ( matrix, 1)
        fprintf(fid,'%10g ', rowTitle(i) );
        for j = 1 : size (matrix, 2 )
%             fprintf(fid,'%12.5f ', matrix(i, j));
            fprintf(fid,'%15g ', matrix(i, j));

        end
        
        fprintf(fid,'\r\n');
    end
end
%% -----------------------------------------------
function write2(fid, vector)
fprintf(fid,'           ');
for j = 1 : size ( vector, 2)
    fprintf(fid,'%15g ', vector(j));
end
fprintf(fid,'\r\n');
end