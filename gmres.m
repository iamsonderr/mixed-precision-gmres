clear; clc; close all
% include file
matrices_loader_from_mat_file


disp('start GMRES.');

matrix = load('./matrix_collection/ns3Da.mat');

% A = matrix('ns3Da').Problem.A;
A = matrix.Problem.A;
b = zeros(size(A,1),1);
x0 = zeros(size(A,1),1);
x0(1) = 1;
restart_m = 100;
tol = 1e-10;% specified accuracy radio 

inner_iteration_count = 0;

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
    if rm <= tol% condition of convergence
        break;
    end
    
    % The specified accuary was not achieved, meet the condition to restarting
    x0 = xm;
end
%% test part
clear; clc; close all

m = matfile('./matrix_collection/494_bus.mat');
m1 = m.Problem
m1.A