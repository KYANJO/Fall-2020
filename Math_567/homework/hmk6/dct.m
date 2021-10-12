function x = dct(x)
%  DCT   Discrete Cosine Transform computed using the fast Fourier Transform.
%     X = dct(x) computes the Discrete Cosine Transform (DCT) of the columns of X.
%
%     X_k = 1/N*(1/2*x_0 + \sum_{j=1}^{N-1} x_j*cos(\pi*j*k/N) + 1/2*(-1)^k*x_N)
%
%     k = 0,...,N, where N is the column length of X.
%
%  See also idct.
[m,~] = size(x);

xe = real( fft( [x; x(m-1:-1:2,:)]/(2*(m-1))));
x = xe(1:m,:);

% k = 0:(m-1);
% D = cos((pi*(k'*k))/(m-1));
% D(:,[1 m]) = 0.5*D(:,[1 m]);
% 
% x = 1/(m-1)*(D*x);

end

