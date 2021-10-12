% clear
clc
% Set up some variables for plotting
LW      =  'LineWidth';
lw      =  1;
clr     =  [221  221  221]/255;
xlbl    =  'Re($\xi$)';
ylbl    =  'Im($\xi$)';
intrptr =  'Interpreter';
ltx     =  'Latex';

%Define the unit circle in the complex plane
N = 1000;
th = linspace(0,2*pi,N);
w = exp(1i*th);
% Solution of characteristic equation in terms of $\xi$
f = @(w) 3*(w.^3-w)./(7*w.^2-2*w+1);
%f = @(w) 12*(w.^3 - w.^2)./(23*w.^2 - 16*w + 5);
% Evaluate $f$ at points on the unit circle and then plot 
% the results:
xi = f(w);
plot(xi,'k-', LW, lw) , hold on
fill(real(xi),imag(xi), clr)
plot([min(real(xi)) max(real(xi))],[0  0], 'b--', LW, lw)
plot([0  0], [min(imag(xi)) max(imag(xi))], 'b--', LW, lw)
xlabel(xlbl, intrptr, ltx), ylabel(ylbl, intrptr, ltx)
xlim([min(real(xi))-0.1 max(real(xi))+0.1])
ylim([min(imag(xi))-0.1 max(imag(xi))+0.1])
daspect([1  1  1]), hold on
f = @(w) 12*(w.^3 - w.^2)./(23*w.^2 - 16*w + 5);
% Evaluate $f$ at points on the unit circle and then plot 
% the results:
xi = f(w);
clr1  =  [100  100  100]/150;
plot(xi,'k-', LW, lw) , hold on
fill(real(xi),imag(xi), clr1)
plot([min(real(xi)) max(real(xi))],[0  0], 'r--', LW, lw)
plot([0  0], [min(imag(xi)) max(imag(xi))], 'r--', LW, lw)
xlabel(xlbl, intrptr, ltx), ylabel(ylbl, intrptr, ltx)
xlim([min(real(xi))-0.2 max(real(xi))+0.3])
ylim([min(imag(xi))-0.3 max(imag(xi))+0.3])
daspect([1  1  1]), hold off