% This function reads VARS input file which includes VARS configuration
% © Saman Razavi
function VARS_inp = read_VARS_inp(filename)

raw_text = read_file(filename); 
trimmed_text = text_trim (raw_text);

VARS_inp.outFldr = trimmed_text{5}; % output folder
VARS_inp.numStars = str2num ( trimmed_text{6} ); % number of center points
VARS_inp.grdSize = str2num ( trimmed_text{7} ); % grid size
VARS_inp.IVARS = str2num (trimmed_text{8} )'; % IVARS specifiers
VARS_inp.mdlFile = trimmed_text{9}; % model file
VARS_inp.mdlFldr = trimmed_text{10}; % model folder
VARS_inp.starFile = trimmed_text{11}; % Star centers location file
VARS_inp.SmpStrtgy = trimmed_text{12}; % sampling strategy
VARS_inp.rndSeed = str2num ( trimmed_text{13} ); % random seed
VARS_inp.btsrpFlg = str2num ( trimmed_text{14} ); % bootstrap flag
VARS_inp.btsrpSize = str2num ( trimmed_text{15} ); % bootstrap size
VARS_inp.btsrpCL = str2num ( trimmed_text{16} ); % bootstrap confidence level
VARS_inp.numGrp = str2num (trimmed_text{17} ); % number of groups
VARS_inp.offlineMode = str2num (trimmed_text{18} ); % offline VARS 
VARS_inp.offlineStage = str2num (trimmed_text{19} ); % offline VARS Stage
VARS_inp.calcFreq = min ( str2num (trimmed_text{20} ), VARS_inp.numStars ); % frequency of VARS calculation and writing
VARS_inp.plotFlg = str2num (trimmed_text{21} ); % plotting flag
VARS_inp.mdlOutLength = str2num (trimmed_text{22} ); % Length of time-ordered model output, "1" for a single-output model
VARS_inp.txtReportFlg = str2num (trimmed_text{23} ); % Text Reporting Flag; enter "0" not to write txt report files, or "1" to write txt report files
end
%% 
function raw_text = read_file(filename)
fileID = fopen(filename);
row_num = 0;
while ~feof(fileID)
    row_num = row_num + 1;
    raw_text{row_num, 1} = fgetl(fileID);
end
fclose (fileID);
end
%% 
function trimmed_text = text_trim (raw_text)
num_lines = size (raw_text, 1);
for i = 1 : num_lines
    comment_index = findstr ( raw_text{ i , 1} , '%' );
    temp = raw_text{ i , 1}(1 : comment_index - 1);
    trimmed_text{i, 1} = strtrim ( temp );
end
end