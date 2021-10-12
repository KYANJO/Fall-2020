clc
%clear

f = @(x) sin((pi*x)/2) + (1/2)*sin(2*pi*x);
g0 = @ (t) 0*t;
g1 = @(t) exp(-(pi^2)*t/4);

uexact = @(x,t) exp(-(pi^2)*t/4).*sin((pi*x)/2) + (1/2)*exp(-(4*pi^2)*t).*sin(2*pi*x);

L = 4:8;

m = 2.^L-1;
N = 2.^L;
alp = 1;

tspan = 1;
Error = zeros(5,1);

for i = 1:5
    
    [u,t,x] = bdf2(f,g0,g1,tspan,alp,N(i),m(i));
    u = u(end,:);
    ue = uexact(x,t(end));
    error = RelNorm(u,ue);
    Error(i) = error;
    
end

loglog(m, Error)
title('Error between approximate and exact solution')

p = polyfit(log(m), log(Error),1)


function L2 = RelNorm(U,Ue)

R = (U - Ue).^2;

L2 = sqrt(sum(R)/sum(Ue.^2));
end
