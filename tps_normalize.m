function  [X, centroid, scale] = tps_normalize(x)
[n, d]=size(x);
centroid=mean(x);
x=x-repmat(centroid,n,1);
scale=sqrt(sum(sum(x.^2,2))/n);
X=x/scale;
end


  
