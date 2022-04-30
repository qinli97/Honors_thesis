function K = tps_kernel(x,y)
[n,d] = size(x);
[m,d] = size(y);
K = zeros(n,m);

R = K;
for i = 1:n
    for j = 1:m
        R(i,j) = norm((x(i,:)-y(j,:)));
    end
    
end

if d ==2
    K = R*log(sqrt(R));  
else
    K = -1*R; %U(r) = -r
    
end

end