% main
clear; clc; close all 

% Ax = b
A = pascal(4);  % A is a pascal matrix.
b = [0 0 0 0]';
x0 = [1 0 0 0]';    % x0 is the initial sovler vector for Ax = b.
m = 4;

% Arnoldi iterative process.
% It generates V which is a orthonormal base /
% for krylov subspace with rows the same as rows of x0 and m+1 columns /
% and H which is a upper hessenberg matrix with m+1 rows and m columns.
[V,H] = Arnoldi(A,b,x0,m);

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
[~,real_m] = size(H);
beta_multiply_e1 = zeros(realm+1,1);beta_multiply_e1(1) = beta;
[T,gm] = Givens( H,beta_multiply_e1 )

% y is the shift from the initial x0.
y1 = inv(T)*bk
%y2 = BackUT( T,bk )
% x = x0 + V * y
realSolution = inv(A)*b
sol1 = x0+V*y1
%sol2 = x0+V(:,1:4)*y2