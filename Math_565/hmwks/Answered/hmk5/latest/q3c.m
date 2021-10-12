clc
clear all

A = 1; B = 0.5;
K = 1-(B/A)^2;
f = @(x) sqrt(1-K*(sin(x)^2));
fprime = @(x) -K*cos(x)*sin(x)*((1-K*sin(x)^2))^(-1/2);

a = 0; b = 1;
Te = fprime(1) - fprime(0);
n = 3:8;
Nt = 2.^n;

Tc = zeros(6,1);
H = zeros(6,1);

Xexact = 0.8866251235367069482;

for j = 1:6
    N = Nt(j);
    h = 1/N;
    H(j) = h;
    
    S = 0;
    for i = 1:N+1
        xi = (i-1)*h;
        
        if ((i == 1) || (i == N+1))
            S = S + f(xi);
        else
            S = S + 2*f(xi);

        end
    end

    S = A*S*h/2 - (h^2/12)*Te;
    
    Tc(j) = S;
end

% error
error = abs(Tc - Xexact);
% Slope
p = polyfit(log(H), log(error),1);

fprintf('The slope of the line is p = %f\n', p(1));
% Plot
loglog(Nt,error)
title('Error against N using Trapezoidal corrected')