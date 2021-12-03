%% main
clear; clc; close all
disp('start GMRES.');

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

A = sprandsym(10,0.7);
b = [0 0 0 0 0 0 0 0 0 0]';
x0 = [1 0 0 0 0 0 0 0 0 0]';
restart_m = 10;
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
    %Resize Rm_bar and gm_bar
    Rm = Rm_bar(1:real_m,1:real_m);
    gm = gm_bar(1:real_m);

    % y is the shift from the initial x0.
    % ym = inv(Rm)*gm;% solve directly
    ym = BackwardUpperTriangular( Rm,gm );% solve backward

    % x = x0 + V * y
    % real_solution = inv(A)*b;
    % xm = x0+Vm*ym;
    xm = x0+Vm*ym;
    xm
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
A = [1,1,1;0,1,1;0,0,1];
b = [4,2,1];
x = zeros(3,1);
[x] = BackwardUpperTriangular(A,b);

