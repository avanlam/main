function [ VARS_inp , VARS_out ] = main_VARS_evolution()

currDir = pwd;
one = read_outputVARS(strcat(currDir,'out_1/'));
ten = read_outputVARS(strcat(currDir,'out_10/'));
hundred = read_outputVARS(strcat(currDir,'out_100/'));


% filename = 'VARS_inp_ten.txt';
% VARS_inp = read_VARS_inp(filename);
% funPath = VARS_inp.mdlFldr;
% factors = read_factorSpace(funPath);
% VARS_inp.factors = factors;
% VARS_inp.plotFlg = 0;
% 
% fprintf('A VARS experiment started: number of neurons = %g, number of stars = %g, minimum h = %g. \n', VARS_inp.numNeurons, VARS_inp.numStars, VARS_inp.grdSize);
% 
% VARS_out = VARS(VARS_inp);

end