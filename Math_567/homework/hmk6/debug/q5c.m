
% clear
clc

f = @(x) exp(-20*(x-pi/2).^2) + (3/2)*exp(-20*(x-4*pi/3).^2);
g = @(x) 0*x;

k = 2*pi/400;
ntime = 400;
m = 200;
c = 1;
h = 2*pi/(m+1);
x = 0:h:2*pi;

[u,v] = waveq(f, g, k, c, ntime,m);

ue = f(x);

waterfall(u(ntime,:))
%waterfall(ue)
