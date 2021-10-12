
close all
clear all

f = @(x) exp(-20*(x-pi/2).^2) + (3/2)*exp(-20*(x-4*pi/3).^2);
g = @(x) 0*x;

m = 200; c = 1;

[u,x,t] = wave1(f,g,m,c);
plot(x,u);

%init condition
un = u(end,:);
uexact = f(x)

%error
error = abs(u - uexact);
errorn = max(error)



function [u,x,t] = wave1(f,g,m,c)

Nt = 400;
h = 2*pi/(m+1);
k = 2*pi/Nt;

x = 0:h:2*pi-h;
j = 0:m; t = j*k;

v=(2*cos((pi*k)/(m+1)))-2;
f=-4*cos(2*x);

%obtaining fcap
fcap=dct(f);

%obtaining ucap
ucap=(h^2)*fcap./v;

%obtaining u
uap=idct(ucap);
u = uap;
end
