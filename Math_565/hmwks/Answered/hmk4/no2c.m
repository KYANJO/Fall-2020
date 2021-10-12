%Code solves the fixed point iteration problem g(x) = 0.1x +1 using
%steffensens method.

clear all; 
close all;

%tolerance
tol = 1e-8;

%intial guess 
x0 = 0;

kmax = 20;

%function g(x,y)
g=@(x) 0.1*x+1;

fprintf('Below is the solution for the fixed point problem;\n');

fprintf('    k          x_k                e_n\n');

[xroot, en] = steffensens(g,x0,tol,kmax)

fprintf('We get convergence in one step since xk =xroot is achived only in one step k=1\n');

fprintf('In this case we require only one iteration to converge to the true solution while \n in 1(e) we require atleast k = 8 iterations depending on the the magnitude of log|eo|.\n');