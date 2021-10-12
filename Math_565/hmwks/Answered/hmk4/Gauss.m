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
fprintf('Gauss-Seidel takes k = %3d\n', k);
x = xkp1;

end