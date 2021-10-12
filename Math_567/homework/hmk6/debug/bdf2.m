function [u,t,x] = bdf2(f,g0,g1,tspan,alp,N,m)

h = 1/(m+1);
dt = tspan/(N+1);

beta = 2*alp*dt/3;
lamb = beta/h^2;

x = 0:h:1;
t = 0:dt:tspan;

U0 = f(x);
U0(1) = g0(t(1));
U0(m+2) = g1(t(1));

[u1,t,x] = cnhteq(f,g0,g1,tspan,alp,N,m);
U1 = u1(2,:);

% solution at t0, t1;
U = zeros(m+2,m+2);
U(1,:) = U0;
U(2,:) = U1;

for k = 3: N+2
    U0(1) = g0(t(k-2));
    U0(m+2) = g1(t(k-2));
    
    U1(1) = g0(t(k-1));
    U1(m+2) = g1(t(k-1));
    
    Un(1) = g0(t(k));
    Un(m+2) = g1(t(k));
    
    % RHS of the system
    d = zeros(m,1);
    
    for j = 2:m
        
        if j == 2
            d(j) = (4/3)*U1(j) - (1/3)*U0(j) +lamb*Un(1);
        
        elseif j == m
            d(j) = (4/3)*U1(j) - (1/3)*U0(j);
            
        else
            d(j) = (4/3)*U1(j) - (1/3)*U0(j) +lamb*Un(m+2);
        end
    end


    % Use colum vectors to construct the matrix of LHS of the system

    a = (1+2*lamb)*ones(m,1);
    b = lamb*ones(m-1,1);

    A = diag(a) + diag(-b,1) + diag(-b,-1);

    % Use of the sparse library to transform A

    A = sparse(A);

    % Solution for the interior points

    U1 = A\d;

    % Append of the boundary solutions to the interior solution

    Un = [Un(1), U1', Un(m+2)];

    U(k,:) = Un;
    U0 = U1;
    U1 = Un;
    
end

u = U;