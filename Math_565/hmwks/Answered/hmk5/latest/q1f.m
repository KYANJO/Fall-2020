clc
clear all

f = @(x) 1./(1+25*x.^2);

N = 10;
x = zeros(1,N+1);
y = zeros(1,N+1);

h = 0.2;

for k = 0:N
    x(k+1) = -1 + k*h;
    y(k+1) = f(x(k+1));
end

m = 301;
h1 = 2/m;
xl = -1:h1:1;

p = zeros(m+1,1);

for j = 1:m+1
    
    p(j) = barycentric(y, x, N, xl(j));

end
%pnew = p(2:m);
fk = f(x);
plot(x, fk)
hold on
plot(xl, p)

function p = barycentric(y, xl,N,x)

S1 = 0;
S2 = 0;
w = zeros(N+1,1);

for j = 1:N+1
    prodt = 1;
    for k = 1:N+1
        if (k ~= j)
            prodt = prodt*(xl(j) - xl(k));
        end
    end
    w(j) = 1/prodt; 
    
    S1 = S1 + (w(j)/(x-xl(j)))*y(j);
    S2 = S2 + (w(j)/(x-xl(j)));
end
p = S1/S2;
end