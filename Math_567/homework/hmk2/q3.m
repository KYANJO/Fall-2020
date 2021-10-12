% Here we calculate solution for boundary value problem using FDM
% by applying DST and iDST via FFT

N = 100; a = 1;  b = -2;
j = [1:N-1]; 
h = pi/N;
x_j = h*j;
f = h^2 * tanh(4*sin(pi*j/N))';
v1 =2*a*cos(pi*j/N) +b;
v=v1';

% DST of vector f via FFT

f_hat = dst(f);
u_hat = f_hat./v;

% Calculating the solution  u

u = idst(u_hat);
plot(x_j, u)
xlabel('x_j')
%axis([0 pi -0.14 0.005])
ylabel('u(x_j)')