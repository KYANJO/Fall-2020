%Conjugate gradients

function u = CG(A, b, tol)

n = size(A,1);

% Intial guesss uo
uo = zeros(n,1);

ro = b - A*uo;
po = ro;


for k = 1:10000
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

fprintf('%3d %12.4e\n', k, norm(rk));

u = uk;

end