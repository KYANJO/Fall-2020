
clear all;
close all;

N = 16;
m = 15;
alp = 1;
tspan = 1;

f = @(x) sin(pi*x/2) + 0.5*sin(2*pi*x);

g0 = @(t) 0;
g1 = @(t) exp((-pi^2*t)/4);

%Numerical solution
[u,t,x] = cnhteq(f,g0,g1,tspan,alp,N,m);

waterfall(x,t,u), view(10,70)
axis([0 1 0 1 0 1] ), xlabel x, ylabel t, zlabel u
title('Numerical Solution');


%Exact solution
u = @(x,t) exp(-pi^2*t'/4).*sin(pi*x'/2) + 1/2*exp(-4*pi^2*t').*sin(2*pi*x');
%u(x,t);
%waterfall(x,t,u(x,t))

%error = abs(u(x,1),u)

