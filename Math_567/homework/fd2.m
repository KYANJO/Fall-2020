m = 39;

ffun = @(x,y) -5*pi^2*sin(pi*x).*cos(2*pi*y);  % Laplacian(u) = f
gfun = @(x,y) sin(pi*x).*cos(2*pi*y);          % u = g on Boundary
uexact = @(x,y) gfun(x,y);   
a = 0;
b = 1;
h = (b-a)/(m+1);
tol = 10^(-8);
maxiter = 1000;
%omega = 2/(1+sin(pi*h));
omega = 1;
tic
[u,x,y] = fd2poissonsor(ffun,gfun,a,b,m, maxiter,omega,tol);

toc

% Plot error
figure, set(gcf,'DefaultAxesFontSize',8,'PaperPosition', [0 0 3.5 3.5]),  
mesh(x,y,u-uexact(x,y)), colormap([0 0 0]),xlabel('x'),ylabel('y'), 
zlabel('Error'), title(strcat('Error, h=',num2str(h))); 