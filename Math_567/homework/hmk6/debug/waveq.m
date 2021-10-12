

function [u,v] = waveq(f, g, k, c, ntime,m)

h = 2*pi/(m+1);
x = 0:h:2*pi;
A = [[0 1] ; [c^2 0]];

u0 = f(x);
v0 = g(x);
q0 = [u0;v0];

lamb1 = (k/6)*(-6*A + 2*k*A^2 - (k^2*A^3)/2);
alpha1 = (k/6)*(k*A^2 - (k^2*A^3)/2 + (k^3*A^4)/4);

lamb = 1 + lamb1;
alpha = alpha1/(12*h);

%alp = 6+2*k+(k^2/2);

%lamb = k*(alp*A + (k-k^2)*A*A)/(6*12*h);

U = zeros(ntime, m+2);
V = zeros(ntime, m+2);
U(1,:) = u0;
V(1,:) = v0;

qn = zeros(2,m+2);

for n = 2:ntime
    
    qn(:,1:2) = q0(:,1:2);
    qn(:,m+1:m+2) = q0(:,m+1:m+2);
    
    for j = 3:m
        
        qn(:,j) = lamb*q0(:,j)+alpha*(q0(:,j-2)-8*q0(:,j-1)+8*q0(:,j+1)-q0(:,j+2));
    end
    
    q0 = qn;
    
    U(n,:) = qn(1,:);
    V(n,:) = qn(2,:);
    %qn = 0;
end
u = U;
v = V;