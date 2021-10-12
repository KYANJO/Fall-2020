k = 150;
u0 = 1;
u1 = 0;
a= 0; b = 1;
l = [7:16];
h = 1./(2.^l);
L2N = zeros(10,1);

for i = 1:10
    
    N = 1/h(i);
    j = [1:N-1];
    %x = zeros(N-1,1);
    x1 = h(i)*j;
    x = x1';
    U = U_app(h(i),k,u0,u1);
    Ue = uexact(x,k);
    %Ue = Ue1';
    l1 = RelL2Norm(U,Ue);
    L2N(i) = l1;
    
end

loglog(h, L2N)

pfit = polyfit(log(h),log(L2N),1);
%p = pfit(2);


function U = U_app(h,k,u0,u1)

N = 1/h;
a = 1+((k*h)^2/12);  b = -2+(5/6)*(k*h)^2;
j = [1:N-1]; 

x = h*j;
f = zeros(N-1,1);
f(j) = h^2;
f(1) = f(1)-a*u0; 
f(N-1) = f(N-1) - a*u1;
v =(2*a*cos(pi*j/N) +b)';

% DST of vector f via FFT

f_hat = dst(f);
u_hat = f_hat./v;

% Calculating the solution  u

U = idst(u_hat);
end

function L2 = RelL2Norm(U,Ue)


R = (U - Ue).^2;

L2 = sqrt(sum(R)/sum(Ue.^2));
end

function ue = uexact(x,k)

ue = (1/k^2)+(1-1/k^2)*cos(k*x)-((1/k^2)+(1-1/k^2)*cos(k))*(1/sin(k))*sin(k*x);

end


