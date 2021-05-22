function [ factors ] = read_factorSpace(name)

file_name = strcat('factorSpace', name, '.txt');

facrSpcFile = strcat(file_name);

raw_text = read_file(facrSpcFile);
[number , name ] = text_trim (raw_text);
[ lb , ub ] = extract_bounds ( number );
factors.lb = lb;
factors.ub = ub;
factors.name = name;
factors.numDim = length(lb);
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
function [number , string ] = text_trim (raw_text)
num_lines = size (raw_text, 1);
for i = 2 : num_lines
    if isempty(raw_text{ i , 1}) == true; break; end
    letters = isletter( raw_text{ i , 1} );
    letter_index = find( letters == 1 , 1 );
    comment_index = findstr ( raw_text{ i , 1} , '%' );

    if isempty(letter_index) == true; letter_index = inf; end
    if isempty(comment_index) == true; comment_index = inf; end
    
    if letter_index == inf && comment_index == inf
        number_untrimmed = raw_text{ i , 1}(1 : end);
        number{i - 1, 1} = strtrim ( number_untrimmed );
        string{i - 1, 1} = 'no name';
    elseif letter_index < comment_index
        number_untrimmed = raw_text{ i , 1}(1 : letter_index - 1);
        number{i - 1, 1} = strtrim ( number_untrimmed );
        if comment_index == inf
            string_untrimmed = raw_text{ i , 1}(letter_index : end);
        else
            string_untrimmed = raw_text{ i , 1}(letter_index : comment_index - 1);
        end
        string{i - 1, 1} = strtrim ( string_untrimmed );
    else
        number_untrimmed = raw_text{ i , 1}(1 : comment_index - 1);
        number{i - 1, 1} = strtrim ( number_untrimmed );
        string{i - 1, 1} = 'no name';
    end
end
end
%%
function [ lb , ub ] = extract_bounds ( number )
num_lines = size(number, 1);
for i = 1 : num_lines
    temp  = str2num (number{i} );
    lb(i, 1) = temp(2); ub(i, 1) = temp(3);
end
end