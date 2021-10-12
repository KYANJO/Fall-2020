clc
clear all

A = 1; B = 0.5;
K = 1-(B/A)^2;
g = @(x) sqrt(1-K*(sin(x)^2));

n = 3:8;
Nt = 2.^n;

Sol = zeros(6,1);
H = zeros(6,1);

Xexact = 0.8866251235367069482;

% Solution for N = 8,16,32,64,128,256
for j = 1:6
    N = Nt(j);
    h = 1/N;
    H(j) = h;
    
    % Trapezoidal rule
    S = 0;
    for i = 1:N+1
        xi = (i-1)*h;
        
        if ((i == 1) || (i == N+1))
            S = S + g(xi);
        else
            S = S + 2*g(xi);

        end
    end

    S = A*S*h/2;
    
    Sol(j) = S;
end

% error
error = abs(Sol - Xexact);

% Slope
p = polyfit(log(H), log(error),1);

fprintf('The slope of the line is p = %f\n', p(1));

% plot

loglog(Nt,error)
title('Error against N using Trapezoidal')

