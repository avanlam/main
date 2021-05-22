function [ indices ] = read_outputVARS(funPath)
filename = 'VARS_out_100.txt';
currDir = pwd;
cd(funPath);
raw_text = read_file(filename);
cd(currDir);

SOBOL_line = 59;
VARS_line = 49;

indices.sobol = text_trim(raw_text{SOBOL_line});
indices.vars = text_trim(raw_text{VARS_line});
indices.vars = indices.vars(2:end);

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
function numbers = text_trim (raw_line)

    trimmed_text = strtrim ( raw_line );
    numbers = str2num(trimmed_text);

end
%%
function [ lb , ub ] = extract_bounds ( number )
num_lines = size(number, 1);
for i = 1 : num_lines
    temp  = str2num (number{i} );
    lb(i, 1) = temp(2); ub(i, 1) = temp(3);
end
end