% Here we calculate solution for boundary value problem using FDM
% by applying DST and iDST via FFT

N = 1000;
u0 = 1; u1 = 0;
h = 1/N;
k = 150;
a = 1+((k*h)^2/12);  b = -2+(5/6)*(k*h)^2;
j = [1:N-1]; 

x = h*j;
f = zeros(N-1,1);
f(j) = h^2;
f(1) = f(1)-a*u0; 
f(N-1) = f(N-1) - a*u1;
v =(2*a*cos(pi*j/N) +b)';

% DST of vector f via FFT

f_hat = dst(f)
u_hat = f_hat./v;

% Calculating the solution  u

u = idst(u_hat);
ue = uexact(x,k);

plot(x, u)
hold on
plot(x, ue)

xlabel('x_j')
ylabel('u(x_j)')

legend('Computed','Exact')

function ue = uexact(x,k)

ue = (1/k^2)+(1-1/k^2)*cos(k*x)-((1/k^2)+(1-1/k^2)*cos(k))*(1/sin(k))*sin(k*x);

end