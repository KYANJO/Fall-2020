% function abm4 solves numerically the general vector-valued IVP.
% It takes in;
% f(t,u) : a vector valued function
% a , b : end points
% uo : initial condition 
% N : number of steps
% abm4 outputs;
% t : vector containing all timesteps
% y : matrix containing the numerical solution of all the components of the
% system u at each timestep.

function [t,u] = abm4(f,a,b,uo,N)

% using RK4 to obtain u1, u2, and u3
k = (b-a)/N;

t = zeros(N+1,1);
u = zeros(N+1,1);

t(1) = a;
u(1) = uo;

for j = 1:3
    u1 = u(j);
    u2 = u(j) + (1/2)*k*f(u1,t(j));
    u3 = u(j) + (1/2)*k*f(u2, (t(j) +k/2));
    u4 = u(j) + k*f(u3, (t(j) +k/2));

    t(j+1) = a + j*k;
    
    u(j+1) = u(j) + (k/6)*(f(u1,t(j)) + 2*f(u2, (t(j) + k/2))...
        + 2*f(u3, (t(j) +k/2)) + f(u4, (t(j) + k)));
end

for n = 4:N
    uast = u(n) + (k/24)*(55*f(t(n),u(n)) - 59*f(t(n-1),u(n-1))...
        + 37*f(t(n-2),u(n-2)) - 9*f(t(n-3),u(n-3)));
    
    t(n+1) = a + n*k;
    
    u(n+1) = u(n) + (k/24)*(9*f(t(n-1),uast) + 19*f(t(n),u(n))...
        - 5*f(t(n-1),u(n-1)) + f(t(n-2),u(n-2)));
end   
     
end    
