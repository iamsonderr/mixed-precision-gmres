%% main
clear; clc; close all
disp('start GMRES.');

matrices_loader_from_mat_file

A = matrix('494.bus').Problem.A;
b = zeros(size(A,1),1);
x0 = zeros(size(A,1),1);
x0(1) = 1;
restart_m = 100;
tol = 1e-6;% specified accuracy radio 

while true

    % Arnoldi iterative process.
    % Input: restart_m is the restart parameter
    [Vm,Hm_bar] = Arnoldi(A,b,x0,restart_m);

    r0 = b-A*x0;
    beta = norm(r0);
    [~,real_m] = size(Hm_bar);
    beta_e1 = zeros(real_m+1,1);beta_e1(1) = beta;
    [Rm_bar,gm_bar] = Givens( Hm_bar,beta_e1 );
    % resize Rm_bar and gm_bar
    Rm = Rm_bar(1:real_m,1:real_m);
    gm = gm_bar(1:real_m);

    % y is the shift from the initial x0.
    % ym = inv(Rm)*gm;% solve directly
    ym = BackwardUpperTriangular( Rm,gm );% solve backward

    % x = x0 + V * y
    % real_solution = inv(A)*b;
    % xm = x0+Vm*ym;
    xm = x0+Vm*ym;
    % judge whether to restart
    rm = norm(b-A*xm);
    if rm/beta <= tol
        break;
    end
    
    % The specified accuary was not achieved, meet the condition to restarting
    x0 = xm;
end
%% test part
clear; clc; close all

file_folder_with_specific_filetype = fullfile('./matrix_collection/','*.mat');% the 1st parameter is file folder name and the 2nd parameter is specific filetype.
structure_of_matching_file_information = dir(file_folder_with_specific_filetype);% read matching file information from choosen file folder with specific filetype.
matching_file_names = {structure_of_matching_file_information.name}% {} is used to generate a cell i.e. yuanzu in chinese.