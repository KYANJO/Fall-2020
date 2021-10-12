
function spline_test_d()

% Function to interpolate

% Data at equispaced points
[filename directory_name] = uigetfile('*.dat', 'Select a file');
fullname = fullfile(directory_name, filename);
A = load(fullname); 

k = A(:,1);
Xdata = A(:,2);
Ydata = A(:,3);

% Build the spline.  Currently, the derivatives at nodes are all set to
% zero.  Your job is to come up with a nicer spline by modifying
% the routines below.
math465 = false;
if (math465)
    pp = math465_build_spline(xdata,Ydata);    
else
    end_cond = 'natural';                         % 'natural' or 'clamped'
    end_data = [Xdata(1); Ydata(end)]';   % For 'clamped' endpoint condition
    pp1 = math565_build_spline(k,Xdata,end_cond,end_data);  % Xdata
    pp2 = math565_build_spline(k,Ydata,end_cond,end_data);  % Ydata
end

% Evaluate the spline at points used for plotting
kv = linspace(1,150,2000);
Xv = ppval(pp1,kv);       % Matlab function

Yv = ppval(pp2,kv);       % Matlab function

% Plot results
figure(4)
clf;

plot(kv,Xv,'r','linewidth',2);
hold on
plot(k,Xdata,'k.','markersize',10);
plot(kv,Yv,'b','linewidth',2);
plot(k,Ydata,'g.','markersize',10);
legend('Spline solution Xdata','XData','Spline solution Ydata', 'Ydata','fontsize',16);

xlabel('x','fontsize',16);
ylabel('y','fontsize',16);
title('Spline interpolation','fontsize',18);

shg

end