% clear
clc
a = 0;
b = 2*pi;
tspan = [0 2*pi];
c = 1;
m = 200;


f = @(x) exp(-20*(x-pi/2).^2) + (3/2)*exp(-20*(x-4*pi/3).^2);
g = @(x) 0*x;
ntime = 400;

[u,v,x] = waveq(f, g, c, ntime,m, tspan, a,b);

ue = f(x);

% Error
u1 = u(end,:);

error1 = abs(u1-ue);
maxerror = max(error1);

% Plot
%waterfall(u(end,:)), view(5,30)
%hold on
%waterfall(ue), view(5,30)
%title('Waterfall plot for two-way wave equation solution')

plot(x, error1)