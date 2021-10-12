
clear all;
close all;

%rng('default')
B = rand(100,100);
A = B'*B;
b = rand(100,1);
tol = 10^(-8);   %relative residual

kmax = 1000;

u = CG(A, b, tol,kmax);

plot(u);
xlabel('dim(u)');
ylabel('u');
title('A graph of solution u');
%Conjugate gradients

function u = CG(A, b, tol,kmax)

n = size(A,1);

% Intial guesss uo
uo = zeros(n,1);

ro = b - A*uo;
po = ro;


for k = 1:kmax
    wo = A*po;
    alphao = (ro'*ro)/(po'*wo);
    uk = uo + alphao*po;
    rk = ro - alphao*wo;
    
    if norm(rk,2)<tol*norm(b,2)
        break;
    end 
    
    betao = (rk'*rk)/(ro'*ro);
    pk = rk + betao*po;
    
    uo = uk;
    ro = rk;
    po = pk;
end

fprintf('The number of iterations k = %3d\n', k);

u = uk;

end