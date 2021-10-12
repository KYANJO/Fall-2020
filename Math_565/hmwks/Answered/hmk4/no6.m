clear all; 
close all;

%rng('default')
B = rand(4,4);
A = B'*B; %To make A symmetrically positive definite
b = rand(4,1);
tol = 10^(-8);   %relative residual
kmax = 10^5;


u = Gauss(A,b,tol,kmax)

%Eigen values of A
lambda = eig(A)

fprintf("Since the eigen vlaues of A are all positive hence A is symmetric positive definite,\n hence given any vector x , x'Ax > 0.\n")

function x = Gauss(A,b,tol,kmax)

n = size(A,1);

% Intial guesss 
xk = zeros(n,1);

%compute an intial residual
rk = b - A*xk;

D = diag(diag(A));
L = tril(A) -D;
M = D + L;

zk = M\rk; %intial approx

for k = 1:kmax
   
    xkp1 = xk +zk;
    rkp1 = b - A*xkp1;
    zkp1 = M\rkp1;
    
    if norm(zkp1) < tol
        break;
    end
    xk = xkp1;
    zk = zkp1;
       
end
fprintf('The number of iterations k = %3d\n', k);
x = xkp1;

end