clear; clc; close all

matrix = load('./matrix_collection/ns3Da.mat');
A = matrix.Problem.A;
b = zeros(size(A,1),1);
x0 = zeros(size(A,1),1);
x0(1) = 1;

restart_m = 100;
tol = 1e-10;% specified accuracy radio 

inner_iteration_counts = 0;
execution_time = 0;

[inner_iteration_counts,execution_time] = GmresUnpreconditionedDouble(A,b,x0,restart_m,tol);
inner_iteration_counts
execution_time
