


function [u,t,x] = cnhteq(f,g0,g1,tspan,alp,N,m)

h = 1/(m+1);
dt = tspan/(N+1);
lamb = (alp*dt)/(2*h^2);

x = 0:h:1;
t = 0:dt:tspan;

U0 = f(x);
U0(1) = g0(t(1));
U0(m+2) = g1(t(1));

% solution at t0;
U = zeros(m+2,m+2);
U(1,:) = U0;

for k = 2: N+2
    U0(1) = g0(t(k-1));
    U0(m+2) = g1(t(k-1));
    Un(1) = g0(t(k));
    Un(m+2) = g1(t(k));
    
    % RHS of the system
    d = zeros(m,1);
    
    for j = 2:m
        
        if j == 2
            d(j) = lamb*U0(j-1) + (1-2*lamb)*U0(j) + lamb*U0(j+1)+lamb*Un(1);
        
        elseif j == m
            d(j) = lamb*U0(j-1) + (1-2*lamb)*U0(j) + lamb*U0(j+1)+lamb*Un(m+1);
            
        else
            d(j) = lamb*U0(j-1) + (1-2*lamb)*U0(j) + lamb*U0(j+1);
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
    U0 = Un;
end

u = U;

        h = 1/(m+1);
dt = tspan/(N+1);
lamb = (alp*dt)/(2*h^2);

x = 0:h:1;
t = 0:dt:tspan;

U0 = f(x);
U0(1) = g0(t(1));
U0(m+2) = g1(t(1));

% solution at t0;
U = zeros(m+2,m+2);
U(1,:) = U0;

for k = 2: N+2
    U0(1) = g0(t(k-1));
    U0(m+2) = g1(t(k-1));
    Un(1) = g0(t(k));
    Un(m+2) = g1(t(k));
    
    % RHS of the system
    d = zeros(m,1);
    
    for j = 2:m
        
        if j == 2
            d(j) = lamb*U0(j-1) + (1-2*lamb)*U0(j) + lamb*U0(j+1)+lamb*Un(1);
        
        elseif j == m
            d(j) = lamb*U0(j-1) + (1-2*lamb)*U0(j) + lamb*U0(j+1)+lamb*Un(m+1);
            
        else
            d(j) = lamb*U0(j-1) + (1-2*lamb)*U0(j) + lamb*U0(j+1);
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
    U0 = Un;
end

u = U;
end

        