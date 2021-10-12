clc
clear all

h = 30;
m = 100;
N = 100;
%a = 0; b = 3000;
t0 = 0.07;
x0 = 1500;
z0 = 1500;
fM = 16;
gamma = 5*10^6;

%vx1 = [-1/2 1/2];
vx1 = [-5/2 -3/2 -1/2 1/2 3/2 5/2];
wx1 = weights(0, vx1,1);

wx2 = wx1(2,:);
wx = [wx2(4:end)  zeros(1,m-6) wx2(1:3)];


Dz = (1/h)*sparse(circulant(wx,1)); % 6 order derivative operator
Dx = -Dz';

% time 
tspan = [0, 0.45];
k = (tspan(2)-tspan(1))/N;
t = 0:k:0.45;
rho = 2200;
% half time step
th = zeros(1,m);

for n = 0:N
    th(n+1) = (n+1/2)*k;
end
% space
x = zeros(m,1);
z = zeros(m,1);
xh = zeros(m,1);
zh = zeros(m,1);
for i = 1:m
    xh(i) = (i-1/2)*h;   % half step for x
    zh(i) = (i-1+1/2)*h; % half step for z
    x(i) = (i-1)*h;
    z(i) = (i-1)*h;
end

for i = 1:size(x,1)
    if x(i) == 1650
        xd = i;
    end
    if z(i) == 1410
        zd = i;
    end
end

% meshgrids 
[X,Z] = meshgrid(x,z);
[Xh,Zh] = meshgrid(xh,zh);   % half step mesh grid
[Xux,Zux] = meshgrid(x,zh);
[Xuz,Zuz] = meshgrid(xh,z);   

lam = zeros(m);     % lame parameter lambda
mn1 = zeros(m);     % lame parameter mu at half space step
mn2 = zeros(m);     % lame parameter mu 
P = zeros(m);

for i = 1:m
lam(i,:) = lambda(Xh(i,:),Zh(i,:));
mn1(i,:) = mu(Xh(i,:),Zh(i,:));
mn2(i,:) = mu(X(i,:),Z(i,:));
P(i,:) = lam(i,:) + 2*mn1(i,:);
end

% initial values
vx0 = zeros(m);
vz0 = zeros(m);
Txx0 = zeros(m);
Tzz0 = zeros(m);
Txz0 = zeros(m);

% computation of ux and uz
ux = zeros(m,m);
uz = zeros(m,m);
uxd = zeros(N,1);
uzd = zeros(N,1);

for n  = 2:N+1
    
    % source +time function
    A = zeros(m);
    for i = 1:m
    
        A(i,:) = Sourcefunction(t(n), X(i,:),Z(i,:),h) ...
                   - Sourcefunction(t(n-1), X(i,:),Z(i,:),h);

    end
    
    % computation of the solutions
    vx0 = vx0 + (k/rho)*(Txx0*Dx' + Dz*Txz0);
    vz0 = vz0 + (k/rho)*(Txz0*Dx' + Dx*Tzz0);
    Txx0 = Txx0 + k*P.*vx0*Dz' + k*lam.*Dz*vz0 + A;
    Tzz0 = Tzz0 + k*lam.*vx0*Dz' + k*P.*Dz*vz0 + A;
    Txz0 = Txz0 + k*mn2.*(vz0*Dz' + Dx*vx0);
    % Txz0 = Txz0 + k*mn2.*(Dz'*vz0 + vx0*Dx);
    % computation of ux and uz
    ux = ux + k*vx0;
    uz = uz + k*vz0;
    
    % ux and uz at the grid (1650, 1410)
    uxd(n) = ux(xd,zd);
    uzd(n) = uz(xd,zd);
    
    % snapshot
    figure(1)
    p = pcolor(Xux,Zux,ux);
    set(p, 'EdgeColor', 'none');
    colormap(gray(100))
    colorbar;
    title(['Time = ',sprintf('%.4f',t(n)),' sec']);
    drawnow;
    
    
end
title('ux')

figure(2)
plot(t, uxd)
hold on
plot(t, uzd)
legend('ux', 'uz')

%lame parameter lambda
function lamb = lambda(x,z)
rho = 2200;
n = size(x,1);
for i = 1:n
cp(i) = Cp(x(i),z(i));
cs(i) = Cs(x(i),z(i));
end

lamb = rho*cp.^2 - 2*rho*cs.^2;
end
% lame parameter mu
function mn = mu(x,z)
rho = 2200;
n = size(x,1);
for i = 1:n
cs(i) = Cs(x(i),z(i));
end

mn = rho*cs.^2;
end

% P wave velocity
function cp = Cp(x,z)

if ((1500 <= x) && (x <= 2100)) && ((1700 <= z) && (z <= 1800))
    
    cp = 1450;
    
else
    
    cp = 3200;
end
end      
% S wave velocity
function cs = Cs(x,z)

if ((1500 <= x) && (x <= 2100)) && ((1700 <= z) && (z <= 1800))
    
    cs = 0;
    
else
    
    cs = 1847.5;
end
end

% source function
function A = Sourcefunction(t, x,z,h)
t0 = 0.07;
x0 = 1500;
z0 = 1500;
fM = 16;
gamma = 5*10^6;

A1 = (1 - 2*pi^2*fM^2*(t-t0)^2);
A2 = exp(-pi^2*fM^2*(t-t0)^2);
Sh1 = Deltafunction(x-x0,h);
Sh2 = Deltafunction(z-z0,h);

A = h^2*gamma*A1*A2*Sh1.*Sh2;

end

% discrete delta function
function Sh = Deltafunction(s,h)

n = size(s,2);

for i = 1:n
    if abs(s(i)) <= 2*h

        Sh(i) = (1/(4*h))*(1+cos((s(i)*pi)/(2*h)));

    elseif abs(s(i)) > 2*h

        Sh(i) = 0;

    end
end
end