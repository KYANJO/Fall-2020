%Using the tridiagonal solver that we developed in problem 3 to solve the
%Helmholtz equation.

k=150;
u0=1;
u1=0;
N=1000;
h=1/N; 
j=[1:N-1]';
a=(1+(1/12)*(k*h)^2)/(h^2);
b=(-2+(5/6)*(k*h)^2)/(h^2);

xj=j/N;
f=zeros(N-1,1);
f(j)=1;

%boundary condition
f(1)=f(1) - a*u0;
f(N-1)=f(N-1) - a*u1;

%obtaining fcap
fcap=dst(f);

%Obtaining ucap
uc=2*a*cos(pi*j/N) + b;
ucap=fcap./uc;

%obtaing u from ucap
u=idst(ucap);
%fprintf('%10s %16.8e\n',u);


%ploting the solution of u and u_ex
plot(xj,u)
hold on
uexac=u_ex(xj,k)
plot(xj,uexac)
legend( 'Numerical','Exact')
ylabel('u(x)')
xlabel('x')
title('A graph of u against x')

%exact
function uexact=u_ex(xj,k)
c=1/k^2;
uexact=c+(1-c)*cos(k*xj)-(c+(1-c)*cos(k))*(csc(k))*sin(k*xj);
end


