function pN = Barycentric(l,y,x,N)

%weights
w = zeros(N+1,1);

sumN = 0; %Numerator
sumD = 0; %Denominator
for j = 1:N+1
    product = 1;
    for k = 1:N+1
        if (k ~= j)
            product = product*(l(j) - l(k));
        end
    end
    w(j) = 1/product; 
    
    sumN = sumN + (w(j)/(x-l(j)))*y(j);
    sumD = sumD + (w(j)/(x-l(j)));
end
pN = sumN/sumD;
end