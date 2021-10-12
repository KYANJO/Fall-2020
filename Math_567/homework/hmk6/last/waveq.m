% This code solve a 1D heat equation
% utt = c^2*uxx
% with initials conditions:
%       at x = 0 : u(x,0) = f(x)
%       at x = 1 : ut(x,0) = g(x)
% using RK4 and fourth-order accurate centered FD
%   Inputs:
%   ntime: number of time step
%   tspan: the time span
%  [a, b]: domain
%   Output:
%     u,v: solutions to the system at each time step
%       x: vector containing spatial discretization

function [u,v,x] = waveq(f, g, c, ntime,m, tspan, a,b)

z = zeros(1,m-3);
v = [0,8,-1,z,1,-8];
D = circulant(v,1);
B = zeros(m+2);

h = (b-a)/(m+1);
k = (tspan(2)-tspan(1))/ntime;
A = (1/(12*h))*[B -D; -c^2*D B];

x = 0:h:2*pi;
H = (k/6)*(6*A+3*k*A^2+k^2*A^3+k^3*A^4/4);

u0 = f(x);
v0 = g(x);
W0 = [u0,v0];
N = size(A,1);
I = eye(N);
G = I+H;
G = sparse(G);

U = zeros(ntime, m+2);
V = zeros(ntime, m+2);

U(1,:) = u0;
V(1,:) = v0;

for n = 2:ntime
    
    W = G*W0';
    u0 = W(1:N/2)';
    v0 = W(N/2+1:end)';
    U(n,:) = u0;
    V(n,:) = v0;
    W0 = [u0, v0];
    
end
u = U;
v = V;
end