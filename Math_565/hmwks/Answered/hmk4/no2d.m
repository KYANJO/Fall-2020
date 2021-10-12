%Code accelerates the convergence of a fixed-point algorithm using
%steffensens method.

clear all; 
close all;

%tolerance
tol = 1e-8;

%intial guess 
x0 = 0.5;

kmax = 100;

%function g(x,y)
g=@(x) (3+3*x-x^2)^(1/3);

fprintf('Below is the solution for the root finding problem;\n');

fprintf('    k          x_k                e_n\n');

[xroot, en] = steffensens(g,x0,tol,kmax)

%Computing e_n
en0 = [];
for k = 1:length(en)-1
   en3 = en(k);
   en0 = [en0,en3];
end
 
%computing e_n+1
en1 = [];
for k = 2:length(en)
   en2 = en(k);
   en1 = [en1,en2];
end
 
figure(1);
loglog(en0,en1);
title("A graph of e_n_+_1 against e_n");
ylabel("e_n_+_1");
xlabel("e_n");


slope_steffensens=polyfit(log(en0),log(en1),1);
slope_steffensens = slope_steffensens(1);
fprintf('slope_steffensens = %f\n',slope_steffensens(1));
fprintf('Hence the steffensens is quadratically convergent since its slope is approximately 2.\n');
%fixed point
[en] = fixed_point(g,x0,tol,kmax);

%computing e_n
enf0 = [];
for k = 1:length(en)-1
   en3 = en(k);
   enf0 = [enf0,en3];
end

%computing e_n+1
enf1 = [];
for k = 2:length(en)
   en2 = en(k);
   enf1 = [enf1,en2];
end

 hold on
 
loglog(enf0,enf1);
legend('Steffensen','fixed point')

slope_fixed_point=polyfit(log(enf0),log(enf1),1);
slope_fixed_point = slope_fixed_point(1);
fprintf('slope_fixed_point = %f\n',slope_fixed_point(1));
fprintf('Hence the fixed point is linearly convergent since its slope is approximately 1.\n');

%fixed point algorithm
function [en]=fixed_point(g,x0,tol,kmax)

xk = x0;
for k = 1:kmax
    xkp1 = g(xk);
    if abs(xkp1 - xk) < tol
        fprintf('Tolerance achieved\n');
        xroot = xkp1;
        break;
    end
    xk = xkp1;
    en(k) = abs(xkp1 - sqrt(3));
end
%fprintf('\n');
%fprintf('Root is %24.16f\n',xkp1);
%fprintf('Number of iterations : %d\n',k);
    
end