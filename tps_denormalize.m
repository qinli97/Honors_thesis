
function x =tps_denormalize(X, centroid, scale)
[m, d]=size(X);
x=X*scale;         % scale
x=x+repmat(centroid,m,1); % move
end