% The code approximates the function f(x) with a 10th polynomial.

clear all
close all

N = 10; h = 2/N; m = 100; 

xx = linspace(-1,1,m)';

%function
f1 = @(x) 1./(1+25*x.^2);

x = zeros(1,N+1);
f = zeros(1,N+1);
for k = 0:N
   x(k+1) = -1 + k*h;
   f(k+1) = f1(x(k+1));
end

PN = [];
for k = 1:m+1
     pN = Barycentric(x,f,xx,N);
     PN = [PN,pN];
end

plot(x,f1(x),'-*'); grid on;
title('f against x');
xlabel('x');ylabel('f');
hold on
plot(xx,PN);
legend('f(x)','P(x)');


function pN = Barycentric(x,f,xx,N)

%weights
w = zeros(N+1,1);

numer = 0; %Numerator
demon = 0; %Denominator
for j = 1:N+1
    temp = 1;
    for k = 1:N+1
        if (k ~= j)
            temp = temp*(x(j) - x(k));
        end
    end
    w(j) = 1/temp; 
    
    xdiff = xx-x(j);
    temp1 = w(j)./(xdiff);
   
    numer = numer + temp1*f(j);
    demon = demon + temp1;
end
pN = numer./demon;
end
