clear; clc; close all

file_folder_name = './matrix_collection/';
specific_filetype = '*.mat';
file_folder_with_specific_filetype = fullfile(file_folder_name,specific_filetype);% generate a cell in which only 1 string exists.
all_structs_of_matching_file_info = dir(file_folder_with_specific_filetype);% read matching file information from choosen file folder with specific filetype.
[row,column] = size(all_structs_of_matching_file_info);

matrix = containers.Map();
for i = 1:row
    matching_matrix_file_name = all_structs_of_matching_file_info(i).name;
    matrix_file_name = strcat(file_folder_name,all_structs_of_matching_file_info(i).name)
    butt_matrix_file_name = matrix_file_name(21:end);
    butt_matrix_file_name_truncation = butt_matrix_file_name(1:end-4);
    matrix(butt_matrix_file_name_truncation) = load(matrix_file_name);
    

end

matrix('apache2').Problem.A

% Ax = b
% A = pascal(4);  % A is a pascal matrix.
% b = [0 0 0 0]';
% x0 = [1 0 0 0]';    % x0 is the initial sovler vector for Ax = b.
% restart_m = 4;

% Ax = b
% A = [1,2,3,4;5,6,7,8;9,10,11,12;13,14,15,16];  % A is a pascal matrix.
% b = [10 26 42 58]';
% x0 = [1 0 0 0]';    % x0 is the initial sovler vector for Ax = b.
% restart_m = 4;

% A = sprandsym(10,0.7);
% b = [0 0 0 0 0 0 0 0 0 0]';
% x0 = [1 0 0 0 0 0 0 0 0 0]';
% restart_m = 10;
% tol = 1e-6;% specified accuracy radio