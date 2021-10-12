% The function numerically solves (2) using fourth-order centered finite differences
% in space and the standard fourth order Runge-Kutta (RK4) method in time.

function [u,v,x,t] = wave(f,g,m,c)

Nt = 400;
h = 2*pi/(m+1);
k = 2*pi/Nt;

x = 0:h:2*pi-h;
j = 0:m; t = j*k;

D1 = (1/(12*h))*(circulant([0,8,-1,zeros(1,m-4),1,-8],1));
A = [zeros(m+1) -D1; -(c^2)*D1 zeros(m+1)];
A = sparse(A);
N = size(A,1);

u0 = f(x);
v0 = g(x);

w0 = [u0,v0];

U = zeros(Nt,m+1); U(1,:) = u0;
V = zeros(Nt,m+1); V(1,:) = v0;

for i = 1:Nt 
    
   wn = w0' + (k/6)*((6*A)+(3*k*A^2)+(k^2*A^3)+(k^3*A^4)/4)*w0';
   u0 = wn(1:N/2)';
   v0 = wn((N/2)+1:end)';
   U(i,:) = u0; V(i,:) = v0;
   w0 = [u0,v0];    
    
end
u=U;
v=V;
end










