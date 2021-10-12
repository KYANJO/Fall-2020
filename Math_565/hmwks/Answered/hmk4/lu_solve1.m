%solves for without pivoting

function x = lu_solve1(L,U,b)

N = size(L,1);

% Forward Solve
y = zeros(N,1);
for i = 1:N
    lk = L(i,1:i-1)';
    yk = y(1:i-1);
    y(i) = b(i) - lk'*yk;
    
end

% Backward Solve
x = zeros(N,1);
for i = N:-1:1
    m = U(i,i);
    xk = x(i+1:end);
    ak = U(i,i+1:end)';
    x(i) = (y(i) - ak'*xk)/m;
end

end