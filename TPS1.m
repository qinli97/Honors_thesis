
%inputs: 

% sample = points from the sample brain
% ref = points from the ref brain

%outputs:
% regist: the transformed sample brain
% w, affine: transformation matrix
% mse: mean standard error of the transformed brain to the ref brain, may
% not need it

function [regist,w,affine, mse] = TPS1(sample, ref)
tic
[n,dim] = size(sample); % n = number of (sample) data points, dim = dimension
[m, dim] = size(ref);   % m = number of (ref) data points, dim = dimension

mse = zeros(n,dim); % construct a zero matrix for mse

[sample, centroid, scale] = tps_normalize(sample); % to normalize the sample brain
[ref, centroid, scale] = tps_normalize(ref); % to normalize the ref brain

% [cell, centroid, scale] = tps_normalize(cell);

K = tps_kernel(sample, sample); % compute the kernel of the sample brain
Pn = [ones(n,1), sample]; % pi = (1, xi, yi, zi), sample brain in homogeneous coordinates

Ls = [K Pn; Pn' zeros(dim+1,dim+1)]; 
param = pinv(Ls)*[ref; zeros(dim+1,dim)]; % transformation parameters

w = param(1:n,:); % wi = (wi, wx, wy,wz), corresponds to different transformation
                  % happens on each axis
affine = param(n+1:end, :); % nonaffine transformations

U = tps_kernel(sample, ref); % distances between sample and ref brains
% [bestlam] = findlam(sample,ref);
% U = U + bestlam*eye(n,n);

P = [ones(m,1), ref]; % ref brain in the homogeneous coordinates
L = [U,P];
regist = L * param; % transformed sample brain


regist = tps_denormalize(regist, centroid, scale); % denormalize the transformed 
                                                    % sample brain to ref
                                                    % brain size
                                                    
ref = tps_denormalize(ref, centroid, scale); % denormalize the ref brain to its 
                                             % original size 

for j= 1:dim
    mse(:,j) = regist(:,j)-ref(:,j); % compute the distances between transformed brain and 
                                     % the ref brain
end
toc

% disp(L)
% disp(param)

figure; % plot the figure
scatter3(regist(:,1),regist(:,2),regist(:,3),'r+');
title('TPS registration');
xlabel('x');
ylabel('y');
zlabel('z');
hold on;
scatter3(ref(:,1),ref(:,2),ref(:,3),'bo');
legend('sample transformed','reference');
% plot3(reconcell(:,1), reconcell(:,2), reconcell(:,3), 'm*', 'MarkerSize',5, 'LineWidth',1);


% label the points
dx = 0.02; dy = 0.02; dz= 0.02;
a = [1:n]'; b = num2str(a); c = cellstr(b);
text(regist(:,1)+dx, regist(:,2)+dy,regist(:,3)+dz,c, 'Color', 'red');
dxx = 0.50; dyy = 0.50; dzz = 0.50;
e = [1:n]'; h = num2str(e); t = cellstr(h);
text(ref(:,1)+dxx, ref(:,2)+dyy,ref(:,3)+dzz,t,'Color','blue');
figure;
% diff = sqrt(sum(regist(:,1)-ref(:,1)).^2+ (regist(:,2)-ref(:,2)).^2+((regist(:,3)-ref(:,3)).^2));
lab = [1:n];
labxyz = {'x','y','z'};

for l = 1:dim
    subplot(3,1,l);

    plot(lab, mse(:,l), 'bo','MarkerSize',10);
    title('Distance from registered sample brain to reference');
    hold on;
    yline(0)
    xlabel('points on the brain');
    ylabel(labxyz{l});
    axis ([1 n -max(max(abs(mse))) max(max(abs(mse)))])
end
%plot(lab, mse, 'bo','MarkerSize',10)
% xlabel('points on the brain');
% ylabel('distance between each correpondence');


% mse = diff; % distance between each pair of correpondence 

end

%%
% function K = tps_kernel(x,y)
% [n,d] = size(x);
% [m,d] = size(y);
% K = zeros(n,n);
% 
% I = eye(n);
% 
% for i = 1:n
%     for j = 1:n
%         K(i,j) = norm((x(i,:)-y(j,:)));
%         K(j,i) = K(i,j);
%     end
%     
% end
% 
% 
% if d ==2
%     K = K*log(sqrt(K));  
% else
%     K = -sqrt(K); %U(r) = -r
%     
% end
% 
% end


% *: mse: diff in x,y,z axis 

%1. lambda: run lambda = 1:100, find the one has lowest mse take log scale

%2. param: 1). fix the param: run the 7 datasets, take average 

%3. fit in the cell points : f(x,y,z) = a0+a1x+a2y+a3z+sum(wi(U(r))); 

%%
% the average of parameters in X2-X6, 9 points

% paramavg = [-0.0050   -0.0016    0.0064
%    -0.0084   -0.0058   -0.0015
%     0.0053    0.0033   -0.0052
%     0.0103    0.0335   -0.0011
%    -0.0042   -0.0033    0.0044
%    -0.0132   -0.0354   -0.0025
%     0.0018    0.0030    0.0031
%     0.0089   -0.0066   -0.0020
%     0.0045    0.0129   -0.0015
%     0.0068    0.0048    0.0017
%     0.8967    0.0866    0.0388
%    -0.0360    1.0429   -0.0993
%    -0.0094   -0.0346    0.7712];



% scaling, centroid, normalize
% distance, for each pooints, in ref brains
% distance for each points in sample brains
% multiply the two to achieve the scaling. 







% fit gaussians to the two data points, uniformly distributed???
% gaussian kernel 

