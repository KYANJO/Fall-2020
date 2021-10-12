
%LU without pivoting
function [L,U] = lu_wopp(A)
    N = size(A,1);
    U = A;
    L = eye(N);     % initializing usingidentity matrix
    
    % Decomposition
    
    for k = 1:N-1
        % Get multiplier, vectors and submatrix
        m = U(k,k);     % multiplier
        ck = U(k+1:end,k);      % column vector 
        ak = U(k,k+1:end)';      % use transpose to get column vector
        
        Ak = U(k+1:end, k+1:end);       % submatrix
        
        % Update L
        
        lk = ck/m;
        L(k+1:end,k) = lk;
        
        % update U
        
        U(k+1:end,k) = 0;       % zero out variables
        U(k+1:end, k+1:end) = Ak - lk*ak';      %  outer product used
        
    end
end