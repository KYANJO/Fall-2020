function pp  = math465_build_spline(xdata,ydata)

N = length(xdata) - 1;
h = diff(xdata);

% TODO : Compute derivatives explicitly.  See page 16-17 of lecture notes.


% This should be replaced with correct derivatives
ddata = zeros(size(xdata));
    
% Get the spline coefficients
coeffs = spline_coeffs(ddata);

% Create spline structure
pp = mkpp(xdata,coeffs);

end


% Spline coefficients for each segment
function coeffs = spline_coeffs(ddata)

global xdata ydata

N = length(xdata) - 1;

h = diff(xdata);
C = [-2, 3, 0, -1; 2, -3, 0, 0; -1, 2, -1, 0; -1, 1, 0, 0];
coeffs = zeros(N,4);
for k = 1:N
    p = [ydata(k); ydata(k+1); ddata(k)*h(k); ddata(k+1)*h(k)];
    S = diag(1./h(k).^(3:-1:0));
    coeffs(k,:) = -p'*C*S;
end

end
