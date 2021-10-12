 %A = [3 -7 -2 2; -3 5 1 0; 6 -4 0 -5; -9 5 -5 12];

 %[L,U,P,pv] = lu_bug_pp(A)

function [L,U,P,pv] = lu_bug_pp(A)

N = size(A,1);

U = A;
L = eye(N);    % Initialize using identity matrix
P = eye(N);
pv = 1:N;

% Decomposition
for k = 1:N-1
    % Find largest pivot in the columnx
    [m,p] = (max(abs(U(k:end,k))));
    U([k,p+k-1],:) = U([p+k-1,k],:);   % Swap rows
    L([k,p+k-1],1:k-1) = L([p+k-1,k],1:k-1);
        
    % Store permutations
    pv([k,p+k-1]) = pv([p+k-1,k]);    
    
    % Get multiplier, vectors and submatrix
     m  = U(k,k);             % Multiplier
    ck = U(k+1:end,k);       % column vector
    ak = U(k,k+1:end)';      % Use transpose to get a column vector
    Ak = U(k+1:end,k+1:end); % Submatrix
    
    % Update L
    lk = ck/m;
    L(k+1:end,k) = lk;   
        
    % Update U
    U(k+1:end,k) = 0;                   % Zero out variables
    U(k+1:end,k+1:end) = Ak - lk*ak';   % Outer product used
end
P = P(pv,:);
end