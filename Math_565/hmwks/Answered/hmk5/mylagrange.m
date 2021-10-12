function Px = mylagrange(xdata,ydata,xvec)
% MYLAGRANGE evaluates interpolating polynomial through points (xdata,ydata)
%
% PX = MYLAGRANGE(XDATA,YDATA,XVEC) where XDATA is a vector of N+1 x-data
% points x0, x1, ... ,x_{N+1} and YDATA is a vector N+1 y-data points y0,
% y1, ... ,y_{N+1}.  The vector XVEC contains x locations at which to
% evaluate the interpolating polynomial
%
% This uses the Lagrange form of the polynomial, which is isn't very
% efficient. A much better way to evaluate the polynomial is to use the
% Barycentric form of the interpolating polynomial.
 

% Return value is the interpolating polynomial evaluated at x
np1 = length(xdata);
n = np1-1;

m = length(xvec);

Pn = zeros(size(xvec));
for k = 1:m  % m = 100 (number of points to plot)
    x = xvec(k);  
    for j = 1:np1    % np1 = 5 
        ljx = lagrange(j,xdata,x);
        Pn(k) = Pn(k) + ljx*ydata(j);
    end
end

end

function lj = lagrange(j,xdata,x)
% Return l_j(x)
%  x     : Scalar value 
%  xdata : N+1 x-data points
%  j     : Index for j'th polynomial

xd = xdata;
xd(j) = [];  % Remove the jth data point. 

num = prod(x - xd);  % (x-xdata(1))*(x - xdata(2))*(x - xdata(3))...

% This part can be computed in part of a 
denom = prod(xdata(j) - xd);

lj = num/denom;


end