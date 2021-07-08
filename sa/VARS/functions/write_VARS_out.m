function write_VARS_out(VARS_inp, VARS_out)
% Author: Saman Razavi in November 2019

% numCntrs   grdSize mdlFile  btsrpSize  btsrpCL
starNum = size( VARS_out.Gamma, 2 );

if exist (VARS_inp.outFldr, 'dir') ~= 7
       mkdir ( VARS_inp.outFldr );
end

if VARS_inp.mdlOutLength == 1 % the model is sinlge output
    filename = strcat (VARS_inp.outFldr, '/VARS_out_', num2str(starNum), '.txt' );
else
    t = VARS_inp.outTime;
    filename = strcat (VARS_inp.outFldr, '/VARS_out_', num2str(starNum),'_t=', num2str(t), '.txt' );
end

fid = fopen(filename,'w');

fprintf(fid,' -------- Variogram Analysis of Response Surfaces (VARS) Output File -------- \r\n');
fprintf(fid,' VARS settings: \r\n');

fprintf(fid,' Target Number of Stars = %g \r\n', VARS_inp.numStars);
fprintf(fid,' Current Number of Stars = %g \r\n', starNum);

fprintf(fid,' minimum h = %g \r\n', VARS_inp.grdSize);

if VARS_inp.btsrpFlg == 1
    fprintf(fid,' Bootstrap Size = %g \r\n', VARS_inp.btsrpSize);
    fprintf(fid,' Bootstrap Confidence Level = %g \r\n', VARS_inp.btsrpCL);
end

hVector = VARS_inp.grdSize * ( 1 :  size ( VARS_out.Gamma{starNum}, 1) )';
% print gamma

fprintf(fid,'\r\n ***  gamma(h) - Directional Variograms: \r\n');
fprintf(fid,'         h        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, hVector, VARS_out.Gamma{starNum});

% print COV
fprintf(fid,'\r\n ***  COV(h) - Directional Covariograms: \r\n');
fprintf(fid,'         h        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, hVector, VARS_out.COV{starNum});

% print ECOV
fprintf(fid,'\r\n ***  ECOV(h) - Directional Expected Covariograms: \r\n');
fprintf(fid,'         h        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, hVector, VARS_out.ECOV{starNum});

% print IVARS
fprintf(fid,'\r\n ***  IVARS(H) - Integrated Variograms: \r\n');
fprintf(fid,'         H        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, VARS_out.IVARSid{starNum}, VARS_out.IVARS{starNum});

% print IVARS based ranks
fprintf(fid,'\r\n ***  Factor Rankings based on IVARS: \r\n');
fprintf(fid,'         H        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, VARS_out.IVARSid{starNum}, VARS_out.rnkIVARS{starNum});

% print VARS-TO
fprintf(fid,'\r\n ***  VARS-TO - Sobol Total-Order Effects calculated through VARS: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
write2(fid, VARS_out.ST{starNum});

% print VARS-TO based ranks
fprintf(fid,'\r\n ***  Factor Rankings based on VARS-TO: \r\n');
fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');

write2(fid, VARS_out.rnkST{starNum});

% print Morris MAEE
fprintf(fid,'\r\n ***  VARS-ABE - Morris mean ABsolute Elementary effects across scales: \r\n');
fprintf(fid,'   delta x        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, hVector, VARS_out.MAEE{starNum});

% print Morris MEE
fprintf(fid,'\r\n ***  VARS-ACE - Morris mean ACtual Elementary effects across scales: \r\n');
fprintf(fid,'   delta x        factor#1        factor#2        factor#3        -------> \r\n');
write1(fid, hVector, VARS_out.MEE{starNum});

if VARS_inp.btsrpFlg == 1
    % print bootstrap results
    fprintf(fid,'\r\n ********         Bootstrapping Results        ******** \r\n');
    % print gamma_LL
    fprintf(fid,' ***  gamma_LL(h) - Lower Limits on Directional Variograms: \r\n');
    fprintf(fid,'         h        factor#1        factor#2        factor#3        -------> \r\n');
    write1(fid, hVector, VARS_out.Gammalb{starNum});
    % print gamma_UL
    fprintf(fid,'\r\n ***  gamma_UL(h) - Upper Limits on Directional Variograms: \r\n');
    fprintf(fid,'         h        factor#1        factor#2        factor#3        -------> \r\n');
    write1(fid, hVector, VARS_out.Gammaub{starNum});
    % print IVARS_LL
    fprintf(fid,'\r\n ***  IVARS_LL(H) - Lower Limits on Integrated Variograms: \r\n');
    fprintf(fid,'         H        factor#1        factor#2        factor#3        -------> \r\n');
    write1(fid, VARS_out.IVARSid{starNum}, VARS_out.IVARSlb{starNum});
    % print IVARS_UL
    fprintf(fid,'\r\n ***  IVARS_UL(H) - Upper Limits on Integrated Variograms: \r\n');
    fprintf(fid,'         H        factor#1        factor#2        factor#3        -------> \r\n');
    write1(fid, VARS_out.IVARSid{starNum}, VARS_out.IVARSub{starNum});
    % print VARS-TO_LL
    fprintf(fid,'\r\n ***  VARS-TO_LL - Lower Limits on Sobol Total-Order Effects calculated through VARS: \r\n');
    fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
    write2(fid, VARS_out.STlb{starNum});
    % print VARS-TO_UL
    fprintf(fid,'\r\n ***  VARS-TO_UL - Lower Limits on Sobol Total-Order Effects calculated through VARS: \r\n');
    fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
    write2(fid, VARS_out.STub{starNum});
    
    % print rlbltyIVARS
    fprintf(fid,'\r\n ***  Reliability Estimates of Factor Rankings based on IVARS: \r\n');
    fprintf(fid,'         H        factor#1        factor#2        factor#3        -------> \r\n');
    write1(fid, VARS_out.IVARSid{starNum}, VARS_out.relIVARS{starNum});
    
    % print VARS-TO based ranks
    fprintf(fid,'\r\n ***  Reliability Estimates of Factor Rankings based on VARS-TO: \r\n');
    fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
    write2(fid, VARS_out.relST{starNum});
    
    % print grouping
    fprintf(fid,'\r\n ***  Factor Grouping based on VARS-TO & IVARS50: \r\n');
    fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
    caption ={'VARS-TO'; 'IVARS50'};
    write3(fid, caption, VARS_out.Gropus{starNum});

    % print rlblty IVARS grps
    fprintf(fid,'\r\n ***  Reliability Estimates of Rankings based on Factor Grouping: \r\n');
    fprintf(fid,'                  factor#1        factor#2        factor#3        -------> \r\n');
    caption ={'VARS-TO'; 'IVARS50'};
    write3(fid, caption, VARS_out.relGrp{starNum});
    
end
fclose(fid);
% end
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
%% -----------------------------------------------
function write3(fid, rowTitle, matrix)
    for i = 1 : size ( matrix, 1)
        fprintf(fid,' %s', rowTitle{i} );
        for j = 1 : size (matrix, 2 )
%             fprintf(fid,'%12.5f ', matrix(i, j));
            fprintf(fid,'%15g ', matrix(i, j));

        end
        
        fprintf(fid,'\r\n');
    end
end