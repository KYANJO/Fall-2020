function x = idct(x)
%  IDCT   Inverse Discrete Cosine Transform computed using the fast Fourier Transform.
%     x = idct(X) computes the inverse Discrete Cosine Transform (DCT) of the columns of X.
%
%     x_j = X_0 + 2\sum_{k=1}^{N-1} X_k*cos(\pi*j*k/N) + 1/2*(-1)^j*X_N,
%
%     j = 0,...,N, where N is the column length of X.
%
%  See also dct.
[m,~] = size(x);
x = (2*m-2)*dct(x);

end

