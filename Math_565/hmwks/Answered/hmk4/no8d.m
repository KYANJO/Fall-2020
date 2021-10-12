%Compares the number of iterations needed using CG and Gauss-Seidel

clear all;
close all;

%rng('default')
B = rand(100,100);
A = B'*B;
b = rand(100,1);
tol = 10^(-8);   %relative residual
kmax = 10^8;


u = CG(A, b, tol,kmax); %Conjugate Gradient
u = Gauss(A,b,tol,kmax); %Gauss-Seidel

fprintf('Therefore CG converges faster than Gauss-Seidel\n');

