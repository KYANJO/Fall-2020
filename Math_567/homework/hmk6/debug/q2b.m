
clc
%clear

f = @(x) sin((pi*x)/2) + (1/2)*sin(2*pi*x);
g0 = @ (t) 0*t;
g1 = @(t) exp(-(pi^2)*t/4);

uexact = @(x,t) exp(-(pi^2)*t/4).*sin((pi*x)/2) + (1/2)*exp(-(4*pi^2)*t).*sin(2*pi*x);

k = 1/16;
h = 1/16;
n = 5;
m = 2^n-1;
N = 2^n-1;
alp = 1;

tspan = 1;

[u,t,x] = cnhteq(f,g0,g1,tspan,alp,N,m);

ue = uexact(x,t);

u = u(m+2,:);
waterfall(u)