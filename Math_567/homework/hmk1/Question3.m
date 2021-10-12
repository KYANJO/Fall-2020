m = input('Enter the number of spaced points, m: ');
% Interval
a = input('Enter the value of a: ');
b = input('Enter the value of b: ');

% space step
h = (b-a)/(m+1);

% Initialization of the matrix of the system
% Mfd is m*m matrix representing the system of the first derivative
% Msd is m*m matrix representing the system of the second derivative
% Ufd and Usd are m*1 column vector repsenting the right hand side of the systems

Mfd = zeros(m,m);
Ufd = zeros(m,1);
Msd = zeros(m,m);
Usd = zeros(m,1);

Mfd(1,1) = 4;
Mfd(1,2) = 1;
Mfd(m,m-1) = 1;
Mfd(m,m) = 4;

Msd(1,1) = 10;
Msd(1,2) = 1;
Msd(m,m-1) = 1;
Msd(m,m) = 10;

for j = 2:1:m-1
    
    xj = a + j*h;
    %first derivative
    
    Mfd(j,j) = 4;
    Mfd(j,j+1) = 1;
    Mfd(j,j-1) = 1;
    
    Ufd(j) = (3/h)*(-u(xj-h) + u(xj+h));
    
    % second derivative
    
    Msd(j,j) = 10;
    Msd(j,j+1) = 1;
    Msd(j,j-1) = 1;
    
    Usd(j) = (12/h^2)*(u(xj-h)-2*u(xj)+u(xj+h));
    
   
end 
Ufd(1) = (3/h)*(-u(a) + u(a+2*h))+up(a);
Ufd(m) = (3/h)*(-u(b-2*h) + u(b))+up(b);

Usd(1) = (12/h^2)*(u(a)-2*u(a+h)+u(a+2*h))+upp(a);
Usd(m) = (12/h^2)*(u(b-2*h)-2*u(b-h)+u(b))+upp(b);

F1 = Mfd\Ufd;
F2 = Msd\Usd;


% Function u we want to approximate its derivative
function u1 = u(x)
    u1 = x^2*exp(-x);
end

% First derivative(u') of the function u
function u1p = up(x)
    u1p = (2*x-x^2)*exp(-x);
end

% Second derivative(u'') of the function u
function u1pp = upp(x)
    u1pp = (x^2-4*x+2)*exp(-x);
end