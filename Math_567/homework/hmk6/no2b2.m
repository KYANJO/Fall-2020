clear all;
close all;

f = @(x) sin(pi*x/2) + 0.5*sin(2*pi*x);

g0 = @(t) 0;
g1 = @(t) exp((-pi^2*t)/4);

%Exact solution
uexact = @(x,t) exp(-(pi^2)*t/4).*sin((pi*x)/2) + 0.5*exp(-(4*pi^2)*t).*sin(2*pi*x);

alp = 1;
tspan = 1;
n = 4:8;
m = 2.^n - 1;
N = 2.^n;

Re_Err = zeros(5,1);

for i =1:5
    
    %Numerical solution
    [u,t,x] = cnhteq(f,g0,g1,tspan,alp,N(i),m(i));
    u = u(end,:);
    %exact solution
    uex = uexact(x,1)';
    Re_Err(i) = RelNorm(u,uex);
    
end

loglog(m,Re_Err);
title('A graph of Relative Error against m');
xlabel('m');
ylabel('Relative error');

p = polyfit(log(m),log(Re_Err),1);
fprintf('The order of accuracy is  %f \n', p(1));

function L2 = RelNorm(U,Uexact)
error = (U - Uexact).^2;
L2 = sqrt(sum(error)/sum(Uexact.^2));
end

