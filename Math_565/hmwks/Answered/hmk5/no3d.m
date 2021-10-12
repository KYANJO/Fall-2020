% This program evaulates the circumference of the ellipse using the trapezoidal rule

clear all
close all

a = 0; b = pi/2; 
A = 1; B = 0.5;
k = sqrt(1 - (B/A)^2);

f = @(x) 4*A*(sqrt(1 - k^2*(sin(x)).^2));


% exact
Tex =  4.84422411027383809921;

Error = [];
C = [];
N = [];
for n = 4:20
    
    N = [N,n];
    Tc = trapezoidal(a,b,f,n);
    C = [C,Tc];
    error = abs(Tc-Tex);
    Error = [Error,error];

end

%log-linear plot
semilogy(N,Error,'-o'); grid on;
title('Error in Circumference calculation vs N');
xlabel('N'); ylabel('Error in Circumference');

% fitting
EN = @(b,N) b(1).*exp(b(2).*N);                                     % Objective Function
B = fminsearch(@(b) norm(Error - EN(b,N)), [-4; -4]) ;                 % Estimate Parameters
figure
plot(N, Error, 'pk'); hold on
plot(N, EN(B,N), '-b')

hold off; grid on

title('Plot of the bestfit line through the data');
xlabel('N'); ylabel('EN and Error in circumference');
text(27, 105, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', B))

%composite trapezoidalrule
function [T] = trapezoidal(a,b,f,n)
    h = (b-a)/n;
    xe = linspace(a,b,n+1); %Nodes at edges
       
    fe = f(xe);
        
    T = (h/2)*(fe(1) + 2*sum(fe(2:end-1)) + fe(end));
end