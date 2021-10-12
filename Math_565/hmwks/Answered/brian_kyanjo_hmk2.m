x = 1;
for k = 1:1075
x = x/2;
if x == 2.163457211e-314
fprintf('%-7d %20.16e\n',k, x);
end 
end
fprintf('\n');
fprintf('%10s %16.8e\n','realmin',realmin);
fprintf('%10s %16.8e\n','x',x);