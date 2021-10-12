% Solution to question 3 (b)
f= @(x,t) (1/pi)*cos(t -x*sin(t));

N=10000; 


x = linspace(0,20,20);
exact = besselj(1,x);
approx = zeros(1,20);
% approximate solutions for N = 8,16,32,64,128,256
  for j=1:20 
    h= pi/N;
    
    t = linspace(0,pi,N+1);
    approx(j) = 0.5*h * (f(x(j),t(1))+2*sum(f(x(j),t(2:N)))+f(x(j),t(end)));
   
  end
% absolute error
error = abs(approx - exact);



% plot

%loglog(N,error)
%title('Error against N using Trapezoidal rule on loglog scale')
%dim = [.3 .5 .3 .3];
%str = "gradient = " +num2str(p(1));
%annotation('textbox',dim,'String',str,'FitBoxToText','on');
%set(gca, 'XDir','reverse')
%grid

% 
plot(x,exact,'r*','markersize',16), hold on
plot(x,approx,'k-')