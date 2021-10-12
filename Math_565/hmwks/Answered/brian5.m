x=logspace(-3,-1,500)
f1=(x.^2)/8
plot(x,f1)
set(gca,’yscale’,’log’)
%f2=abs(((sqrt(1+(x.^2))-1)/(x.^2))-(1/2))
set(gca,'yscale