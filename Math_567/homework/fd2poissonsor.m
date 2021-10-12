% Numerical approximation to Poisson's equation over the square [a,b]x[a,b] with
% Dirichlet boundary conditions.  Uses a uniform mesh with (n+2)x(n+2) total
% points (i.e, n interior grid points).
% Input:
%     ffun : the RHS of poisson equation (i.e. the Laplacian of u).
%     gfun : the boundary function representing the Dirichlet B.C.
%      a,b : the interval defining the square
%        m : m+2 is the number of points in either direction of the mesh.
% Ouput:
%        u : the numerical solution of Poisson equation at the mesh points.
%      x,y : the uniform mesh.
%
function [u,x,y] = fd2poissonsor(ffun,gfun,a,b,m, maxiter,omega,tol)

h = (b-a)/(m+1);   % Mesh spacing

[x,y] = meshgrid(a:h:b);   % Uniform mesh, including boundary points.

idx = 2:m+1;
idy = 2:m+1;
u = zeros(m,m);

% Compute boundary terms, south, north, east, west
ubs = feval(gfun,x(1,1:m+2),y(1,1:m+2));     % Include corners
ubn = feval(gfun,x(m+2,1:m+2),y(m+2,1:m+2)); % Include corners
ube = feval(gfun,x(idy,m+2),y(idy,m+2));     % No corners
ubw = feval(gfun,x(idy,1),y(idy,1));         % No corners

% Evaluate the RHS of Poisson's equation at the interior points.
f = feval(ffun,x(idy,idx),y(idy,idx));

% Adjust f for boundary terms
f(:,1) = f(:,1) - ubw/h^2;             % West
f(:,m) = f(:,m) - ube/h^2;             % East
f(1,1:m) = f(1,1:m) - ubs(idx)/h^2;    % South
f(m,1:m) = f(m,1:m) - ubn(idx)/h^2;    % North

%f = reshape(f,m*m,1);



for k = 0:maxiter
    for j = 2:(m-1)
        for i = 2:(m-1)
            u(i,j) = (1-omega)*u(i,j)+(omega/4)*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)-h^2*f(i,j));
        end
    end
    
    residual = zeros(m,m);
    for j = 2:(m-1)
        for j = 2:(m-1)
            residual(i,j) = -4*u(i,j)+(1-omega)*u(i,j)+(omega/4)*(u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)-h^2*f(i,j));
        end
    end
    
	if norm(residual(:),2)<tol*norm(f(:),2)
		break
	end
end


% Convert u from a column vector to a matrix to make it easier to work with
% for plotting.
%u = reshape(u,m,m);

% Append on to u the boundary values from the Dirichlet condition.
u = [ubs;[ubw,u,ube];ubn];
 
end



