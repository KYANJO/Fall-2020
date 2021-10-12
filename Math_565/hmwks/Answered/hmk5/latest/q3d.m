clc
clear all

A = 1; B = 0.5;
K = 1-(B/A)^2;
g = @(x) sqrt(1-K*(sin(x)^2));
a = 0; b = pi/2;
%n = 4:20;
Nt = 4:20;
l = size(Nt,2);
Sol = zeros(l,1);
H = zeros(l,1);

Xexact = 4.84422411027383809921;

% Solution for N
for j = 1:l
    N = Nt(j);
    h = (b-a)/N;
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

    S = 4*A*S*h/2;
    
    Sol(j) = S;
end

% error
error = abs(Sol - Xexact);

% plot

semilogy(Nt,error)
title('Error against N using Trapezoidal')

