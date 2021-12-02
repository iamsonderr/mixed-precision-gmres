% main
clear; clc; close all 
disp('start GMRES.');
% Ax = b
A = pascal(4);  % A is a pascal matrix.
b = [0 0 0 0]';
x0 = [1 0 0 0]';    % x0 is the initial sovler vector for Ax = b.
restart_m = 4;

% Arnoldi iterative process.
% Input: restart_m is the restart parameter
[Vm,Hm_bar] = Arnoldi(A,b,x0,restart_m);

% Givens rotation to hessenberg matrix and /
% beta * e1 from Arnoldi process, making /
% hessenberg matrix into upper triangular matrix.
% Input: be is beta * e1
%        H is hessenberg matrix
% Output: T is the matrix after Givens rotation /
%               hessenberg matrix
%         bk is the vector after Givens ritation /
%               to beta * e1
r0 = b-A*x0;
beta = norm(r0);
[~,real_m] = size(Hm_bar);
beta_e1 = zeros(real_m+1,1);beta_e1(1) = beta;
[Rm_bar,gm_bar] = Givens( Hm_bar,beta_e1 );
%Resize Rm_bar and gm_bar
Rm = Rm_bar(1:real_m,1:real_m);
gm = gm_bar(1:real_m);

% y is the shift from the initial x0.
y1 = inv(Rm)*gm;% solve directly
y2 = BackwardUpperTriangular( Rm,gm );% solve backward

% x = x0 + V * y
realSolution = inv(A)*b;
solution1 = x0+Vm*y1;
solution2 = x0+V(:,1:4)*y2