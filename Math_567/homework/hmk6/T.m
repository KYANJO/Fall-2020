clear all;
close all;

N = 10;
m = 5;

alp = 1;
k = 1/16;
h = 1/16;
tspan = 1;
%[x,t] = meshgrid(0:h:1);
%t = 0:k:1;

f = @(x) sin(pi*x/2) + (1/2)*sin(2*pi*x);

g0 = @(x) 0*x;

g1 = @(t) exp(-pi^2*t/4);

%Numerical solution
[u,t,x] = cnhteq(f,g0,g1,tspan,alp,N,m);
waterfall(x,t,u)



% This function Numerically solves the 1-D heat equation using the
% Crank-Nicolson scheme
% input
% f : function representing intial condition
% g0 and g1 : function represents boundary conditions
% tspan : time span over
% N : Number of time steps
% m : Number of points for the uniform spatial discretisation.
% Return
% u : Approximate solution at each time step (matrix)
% t : Vector containing all time steps
% x : Vector containing the spatial discretization.

function [u,t,x] = cnhteq(f,g0,g1,tspan,alp,N,m)

h = 1/(m+1);
dt = tspan/(N+1);
r = (alp*dt)/(2*h^2);

x = 0:h:1;
t = 0:dt:tspan;

U0 = f(x);
U0(1) = g0(t(1));
U0(m+2) = g1(t(1));

% solution at t0;
U = zeros(m+2,m+2);
U(1,:) = U0;

for k = 2: N+2
    U0(1) = g0(t(k-1));
    U0(m+2) = g1(t(k-1));
    Un(1) = g0(t(k));
    Un(m+2) = g1(t(k));
    
    % RHS of the system
    d = zeros(m,1);
    
    for j = 1:m
        
        if j == 2
            d(j) = r*U0(j) + (1-2*r)*U0(j+1) + r*U0(j+2)+r*Un(1);
        
        elseif j == m
            d(j) = r*U0(j) + (1-2*r)*U0(j+1) + r*U0(j+2)+r*Un(m+1);
            
        else
            d(j) = r*U0(j) + (1-2*r)*U0(j+1) + r*U0(j+2);
        end
    end



    % Use colum vectors to construct the matrix of LHS of the system

    a = (1+2*r);
    b = r;

    % Use of the sparse library to transform A

    A = sparse(toeplitz([a b zeros(1, m-2)])); 
    
    % Solution for the interior points

    U1 = A\d;

    % Append of the boundary solutions to the interior solution

    Un = [Un(1), U1', Un(m+2)];

    U(k,:) = Un;
    U0 = Un;
end

u = U;

end







