% This code solve a 1D heat equation
% ut = alpha*uxx
% with boundary conditions:
%       at x = 0 : u(0,t) = g0(t)
%       at x = 1 : u(1,t) = g1(t)
% and initial condition:
%       u(x,0) = f(x)
% using BDF2 method.
% Output:
%   u: solution at each time step
%   t: vector containing the time steps
%   x: vector containing spatial discretization

function [u,t,x] = bdf2(f,g0,g1,tspan,alp,N,m)

h = 1/(m+1);
dt = tspan/N;
beta = 2*alp*dt/3;
lamb = beta/h^2;

x = 0:h:1;
t = 0:dt:tspan;

U0 = f(x);

u1 = cnhteq(f,g0,g1,tspan,alp,N,m);
U1 = u1(2,:);

% solution at t0, t1;
U = zeros(m+2,m+2);
U(1,:) = U0;
U(2,:) = U1;
gleft = zeros(m,1);

% Use colum vectors to construct the matrix of LHS of the system

a = (1+2*lamb)*ones(m,1);
b = lamb*ones(m-1,1);

A = diag(a) + diag(-b,1) + diag(-b,-1);

% Use of the sparse library to transform A

A = sparse(A);


for k = 3: N+1
    
    gleft(1) = g0(t(k));
    gleft(m) = g1(t(k));
    
    % RHS of the system
    d = (4/3)*U1(2:m+1)' - (1/3)*U0(2:m+1)' + lamb*gleft;
    
    % Solution for the interior points

    Un = A\d;

    % Append of the boundary solutions to the interior solution

    Un = [gleft(1), Un', gleft(m)];

    U(k,:) = Un;
    U0 = U1;
    U1 = Un;
    
end
u = U;
end
