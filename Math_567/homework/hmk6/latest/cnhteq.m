% This code solve a 1D heat equation
% ut = alpha*uxx
% with boundary conditions:
%       at x = 0 : u(0,t) = g0(t)
%       at x = 1 : u(1,t) = g1(t)
% and initial condition:
%       u(x,0) = f(x)
% using Crank-Nicolson method.
% Output:
%   u: solution at each time step
%   t: vector containing the time steps
%   x: vector containing spatial discretization

function [u,t,x] = cnhteq(f,g0,g1,tspan,alp,N,m)

h = 1/(m+1);    % space step
k = tspan/N;    % time step

U0 = zeros(1,m+2);
x  = zeros(m+2,1);
t  = zeros(m+2,1);

for i = 1:m+2
    x(i) = (i-1)*h;
    U0(i) = f(x(i));  % initial solution
end

for j = 1:N+1
    t(j) = (j-1)*k;
    
end

lamb = (alp*k)/(2*h^2);

% solution at t = 0;
U = zeros(N+1,m+2);
U(1,:) = U0;

gleft = zeros(m,1);
gright = zeros(m,1);

ar = (1+2*lamb)*ones(m,1);
al = (1-2*lamb)*ones(m,1);   
b = lamb*ones(m-1,1);
% The matrix of LHS and RHS of the system
A = diag(ar) + diag(-b,1) + diag(-b,-1); % LHS
B = diag(al) + diag(b,1) + diag(b,-1);   % RHS
% Use of the sparse library to transform A
A = sparse(A);
B = sparse(B);

for n = 1: N
    
    % Boundary conditions
    
    gleft(1) = g0(t(n+1));
    gleft(m) = g1(t(n+1));
    gright(1) = g0(t(n));
    gright(m) = g1(t(n));

    % Solution for the interior points
    d = B*U0(2:m+1)' + lamb*(gleft+gright);
    U1 = A\d;

    % Append of the boundary solutions to the interior solution

    Un = [gleft(1), U1', gleft(m)];

    U(n+1,:) = Un;
    U0 = Un;
end

u = U;
end
