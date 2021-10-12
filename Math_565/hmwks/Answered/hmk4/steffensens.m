

function [xroot, en] = steffensens(g,x0,tol,kmax)


xkm1 = x0;

for k = 1:kmax
   gk = g(xkm1);
   ggk = g(gk);
   D = (ggk - 2*gk +xkm1);
   if (D==0)
       fprintf('Tolerance achived\n')
       xroot = g(xkm1);
       break;
   else
       xk = xkm1 - (gk-xkm1)^2/D;
       
   end
   
   
   en(k) = abs(xk - xkm1);
  
   fprintf('%5d %20.16e, %12.4e\n',k,xk,en(k));
   if (en(k) < tol)
       fprintf('Tolerance achieved\n');
       xroot = xk;
       break;
   end
   xkm1 = xk;
    
end
xroot = xk;
end