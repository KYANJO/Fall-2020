% Solution to question 3 (d)

A = 1; B = 0.5;
K = 1-(B/A)^2;
f = @(x) 4*A*sqrt(1-K*(sin(x).^2));


N = (4:20)';
approx = zeros(17,1);
H = zeros(17,1);

exact = 4.84422411027383809921;
% approximate solutions for N = 4 to 20
for j = 1:17
    
    H(j)= (pi/2)/N(j);
    
    t = linspace(0,pi/2,N(j)+1);
    approx(j) = 0.5*H(j) * (f(t(1))+2* sum(f(t(2:N(j))))+f(t(end)));
    
end

% absolute error
error = abs(approx - exact);

% plot

figure(1)
semilogy(N,error)
title('Error against N using Trapezoidal rule on semilog scale')
set(gca, 'XDir','reverse')
xlabel('N')
ylabel('Log of absolute error')
grid

% fitting line of bestfit to the error
f = @(b,N) b(1).*exp(b(2).*N);                                     % Objective Function
B = fminsearch(@(b) norm(error - f(b,N)), [-4; -4]) ;                 % Estimate Parameters
figure
plot(N, error, 'pg')
hold on
plot(N, f(B,N), '-r')
hold off
grid
xlabel('N')
ylabel('error')
text(27, 105, sprintf('f(x) = %.1f\\cdote^{%.3f\\cdotx}%+.1f', B))