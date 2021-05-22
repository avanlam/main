function [ pivots , X_star , X_scaled_mat, Y_offline ] = read_STAR_samples( VARS_inp )
% for off-line mode
currDir = pwd;
cd ( VARS_inp.outFldr );
load STAR_out;

if VARS_inp.mdlOutLength == 1 % the model is sinlge output
    filename = 'STAR_in.smp';
else
    t = VARS_inp.outTime;
    filename = strcat('STAR_in_t=' , num2str(t), '.smp' );
end
if exist(filename, 'file') == 0
    error('Error: The actual length of model output sequence is shorter than that specifed in VARS_inp.txt.');
end
Y_offline = importdata( filename ); % the result of model runs externally

cd(currDir);
end