%computing u using DST
N=100;
a=1;
b=-2;
h=pi/N;

j=[1:N-1];
xj=h*j;
f=(h^2)*tanh(4*sin(pi*j/N));
ft=transpose(f);
%Obtaining fcap
fcap=dst(ft);

%Obtaining ucap
uc=2*a*cos(pi*j/N) + b;
ucap=fcap./transpose(uc);

%obtaing u from ucap
u=idst(ucap);
fprintf('%10s %16.8e\n',u);

%ploting the solution of u
plot(xj,u,'*');
ylabel('u(x)');
xlabel('x');
title('A graph of u against x');

