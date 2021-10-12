% This program corrects the trapezoidal method to obtain fourth order
% method.

clear all
close all

a = 0; b = 1; 
A = 1; B = 0.5;
k = sqrt(1 - (B/A)^2);

f = @(x) A*(sqrt(1 - k^2*(sin(x)).^2));
fp = @(x) -A*k^2*sin(2*x).*(2*sqrt(1 - k^2*(sin(x)).^2)).^-1; %fprime

%exact solution 
Tex =  0.8866251235367069482;

n = [8,16,32,64,128,256];
c = length(n);
Error = [];

for i = 1:c
    
    Tc= trape(a,b,f,fp,n(i));
    error = abs(Tc-Tex);
    Error = [Error,error];

end

%loglog plot
loglog(n,Error,'-*');
title('Errors vs N');
xlabel('N'); ylabel('Errors');

%order of convergence
p = polyfit(log(n),log(Error),1); p(1)

fprintf('Hence order of convergence is 4\n');


function [Tc] = trape(a,b,f,fp,n)
    h = (b-a)/n;
    xe = linspace(a,b,n+1); %Nodes at edges
       
    fe = f(xe);
        
    T = (h/2)*(fe(1) + 2*sum(fe(2:end-1)) + fe(end));
    
    fpp = fp(xe);
    
    Tc = T - ((h^2)/12)*(fpp(end) - fpp(1));
end