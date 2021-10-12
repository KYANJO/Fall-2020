close all
clear all

f = @(x) exp(-20*(x-pi/2).^2) + (3/2)*exp(-20*(x-4*pi/3).^2);
g = @(x) 0*x;

m = 200; c = 1;

[u,v,x,t] = wave(f,g,m,c);

figure(1);
plot(x,u);
title('Solution on equispaced grid');
xlabel('x'); ylabel('u');

figure(2);
plot(t,u);
title('u against t');
xlabel('t'); ylabel('u');

figure(3)
waterfall(u), view(10,70)
%axis([0 1 0 1 0 1] ), xlabel x, ylabel t, zlabel u
title('Waterfall plot');



%init condition
un = u(end,:);
uexact = f(x);

%error
error = abs(un - uexact);
errorn = max(error)

figure(4);
plot(x,error);
title('Error at t=2*pi');
xlabel('x'); ylabel('Error');
