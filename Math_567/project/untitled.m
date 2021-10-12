% This program uses 6th order approximation of spartial derivatives. It
% numerically solves the elastic wave equation.

clear all
close all

m = 100; N = 100;  %grid cells

%density over the whole domain
rho = 2200; %kg/m^3

%space
a = 0; b = 3000; 
%h = (b-a)/m;
h = 30;
xn = zeros(m,1);
x = zeros(m,1);
zn = zeros(m,1);
z = zeros(m,1);
for i = 1:m
   xn(i) = (i-1)*h;
   x(i) =  (i-0.5)*h; 
   zn(i) = (i-1)*h;
   z(i) =  (i-0.5)*h; 
end

%for i = 1:size(x,1)
%    if x(i) == 1650
%        xd = i;
%    end
%    if z(i) == 1410
%        zd = i;
%    end
%end


%meshes
[X,Z] = meshgrid(x,z);
[Xn,Zn] = meshgrid(xn,zn);   
[Xnx,Znx] = meshgrid(x,zn);
[Xnz,Znz] = meshgrid(xn,z);   

lamb = zeros(m);
mu1  = zeros(m);
mu2 = zeros(m);
P = zeros(m);
for i = 1:m
    lamb(i,:) = lambda(Xn(i,:),Zn(i,:));
    mu1(i,:) = mu(Xn(i,:),Zn(i,:));
    mu2(i,:) = mu(X(i,:),Z(i,:));
    P(i,:) = lamb(i,:) + 2*mu1(i,:); 
end

% reshape 
lamb = reshape(lamb,m*m,1);
mu1 = reshape(mu1,m*m,1);
mu2 = reshape(mu2,m*m,1);
P = reshape(P,(m)*(m),1);

%Dx and Dz
D6 = (1/h)*circulant([0,1.1719,-0.0651,0.0417,zeros(1,m-7),-0.0417,0.0651,-1.1719],1); %derivative operator 6th order
D6 = sparse(D6);
Dx = kron(speye(m,m),D6); Dz = kron(D6,speye(m,m));

% Displacement initial values
vx0 = zeros(m);
vz0 = zeros(m);

% Stress intial Values
Txx0 = zeros(m);
Tzz0 = zeros(m);
Txz0 = zeros(m);   

% at (x,z) = (1650m, 1410m)
%xclose = 1650; zclose = 1410;

%time
at =0; bt = 0.45; 
k = (bt-at)/N; %timestep
t = at:k:bt;
%half
tn = zeros(1,m);
for i = 0:N
   tn(i+1) = (i + 0.5)*k;
end

% ux and uz
ux = zeros(m,m);
uz = zeros(m,m);
uxd = zeros(N,1);
uzd = zeros(N,1);

for n  = 2:N+1  
   
    vx0 = reshape(vx0,(m)*(m),1);
    vz0 = reshape(vz0,(m)*(m),1);
    Txx0 = reshape(Txx0,(m)*(m),1);
    Tzz0 = reshape(Tzz0,(m)*(m),1);
    Txz0 = reshape(Txz0,(m)*(m),1);
    
    A = zeros(m);
    for i = 1:m
        A(i,:) = Sourcetime(t(n), Xn(i,:),Zn(i,:),h);
    end
    
    A = reshape(A,(m)*(m),1);
    
    % computation of the solutions
    vxn = vx0 + (k/rho)*(Dx*Txx0 + Dz*Txz0);
    vzn = vz0 + (k/rho)*(Dx*Txz0 + Dz*Tzz0);
    Txxn = Txx0 + k*P.*Dx*vxn + k*lamb.*Dz*vzn + k*A;
    Tzzn = Tzz0 + k*lamb.*Dx*vxn + k*P.*Dz*vzn + k*A;
    Txzn = Txz0 + k*mu2.*(Dx*vzn + Dz*vxn);
    
    % reshape into matrix form 
    vx0 = reshape(vxn,m,m);
    vz0 = reshape(vzn,m,m);
    Txx0 = reshape(Txxn,m,m);
    Tzz0 = reshape(Tzzn,m,m);
    Txz0 = reshape(Txzn,m,m);
    
    % computation of ux and uz
    ux = ux + k*vx0;
    uz = uz + k*vz0;
    
    %uxd(n) = ux(xd,zd);
    %uzd(n) = uz(xd,zd);
    
    figure(1)
    p = pcolor(Xnx,Znx,ux);
    set(p, 'EdgeColor', 'none');
    xlabel('time(sec)'); ylabel('displacement (millimeters');
    colormap(gray(100))
    colorbar;
    drawnow;
    
end

%time-series plot of ux and uz
%figure(2)
%plot(t', uxd)
%xlabel('time(sec)'); ylabel('displacement (millimeters)');
%tile('Tim-series of recorded displacements near (x,z) = (1650,1410)');
%legend(u_x,u_z);
%hold on
%plot(t', uzd)

%lambda
function la = lambda(x,z)

rho = 2200;
n = size(x,1);

for i = 1:n
Cp(i) = CP(x(i),z(i));
Cs(i) = CS(x(i),z(i));
end

la = rho*Cp.^2 - 2*rho*Cs.^2;
end

% mu
function mhu = mu(x,z)

rho = 2200;
n = size(x,1);

for i = 1:n
Cs(i) = CS(x(i),z(i));
end

mhu = rho*Cs.^2;
end

% Discrete delta function
function delta = Deltah(ep,h)

n = size(ep,2);

for i = 1:n
    
    if abs(ep(i)) <= 2*h
        delta(i) = (1/(4*h))*(1+cos((ep(i)*pi)/(2*h)));
    elseif abs(ep(i)) > 2*h
        delta(i) = 0;
    end
    
end
end

% derivative of source-time function
function S = Sourcetime(t, x,z,h)
t0 = 0.07; %sec
x0 = 1500; %m
z0 = 1500; %m
fM = 16;   %Hz
gamma = 5*10^6; %Pa

S1 = Deltah(x-x0,h);
S2 = Deltah(z-z0,h);

S = h^2*gamma*(-6*pi^2*fM^2*(t-t0)+4*pi^4*fM^4*(t-t0)^3)*(exp(-(pi^2)*(fM^2)*(t-t0)^2))*S1.*S2;

end

%pressure-wave velocity
function cp = CP(x,z)

if ((1500 <= x) && (x <= 2100)) && ((1700 <= z) && (z <= 1800))
        cp = 1450; %m/s
else
        cp = 3200; %m/s
end

end    

%shear-wave velocity
function cs = CS(x,z)

if ((1500 <= x) && (x <= 2100)) && ((1700 <= z) && (z <= 1800))
      cs = 0;   %m/s
else
      cs = 1847.5; %m/s
end

end