clear all;
close all;

J =@(x,t) 1/pi*cos(t -x*sin(t));

N = 20; h = 20/N;  a=0; b=pi; n = 1e6; h1 = (b-a)/n;
t = linspace(0,pi,n+1);

x = [];
T = [];
for i = 0:19
    x1 = i*h;
    x = [x,x1];
end

for i = 1:20
     %trapezoidal rule 
    T1 = (h1/2)*(J(x(i),t(1))+2*sum(J(x(i),t(2:n)))+J(x(i),t(end)));
    T = [T,T1];
end

%exact
J1=besselj(1,x);
error = abs(J1 - T); 

%table
N = [1:20]';
Table = table(N(:),error(:),'VariableNames',{'N','Error'})  

fprintf('Hence error values are approximately 10^-16');

%convergence
fprintf('Exponential convergence, due to the fact that we are dealing with a periodic integral');

figure(1)
loglog(N,error, '-o'); grid on
xlabel('N');ylabel('Error');
title('Error \approx 10^-^1^6 vs N');

figure(2)
plot(x,T); grid on
xlabel('x');ylabel('J');
title('Bessel function vs x');
hold on
plot(x,J1,'-o');
legend('J_t_r_a_p_e_z_i_o_d_a_l','J_e_x_a_c_t')
