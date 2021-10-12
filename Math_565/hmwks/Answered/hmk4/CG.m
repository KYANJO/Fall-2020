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

fprintf('CG takes k = %3d\n', k);

u = uk;

end