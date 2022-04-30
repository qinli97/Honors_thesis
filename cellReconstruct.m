close all;
[regist, w,affine, mse] = TPS1(sample2,ref);
cell = [205.433	-12892.132	2050.563 % single cell coordinates
210.501	-12908.86	2051.716
204.55	-12917.63	2060.658
132.899	-12947.208	2070.649
70.819	-12999.958	2064.005
31.479	-13062.646	2083.867
-29.891	-13131.884	2083.867
-107.518	-13219.57	2083.867
-137.841	-13288.058	2087.197
-150.044	-13297.75	2071.898
-173.989	-13311.45	2050.251
-132.014	-13310.208	2070.896
-142.261	-13365.768	2078.702
-89.341	-13384.852	2104.034
214.118	-12919.568	2049.94
182.4	-12950.886	2023.334
121.186	-12981.456	1992.957
53.573	-13007.99	1967.481
-46.014	-13033.2	1934.529
-138.767	-13067.604	1916.739
-202.569	-13092.92	1902.071
-269.409	-13111.88	1898.638
-334.954	-13111.876	1898.638
-389.348	-13124.938	1861.061
-462.113	-13133.942	1829.057
-491.8	-13125.52	1825.717
225.284	-12916.6	2050.045
265.709	-12900.928	2074.814];

cell(:,2) = -1*cell(:,2); % make sure the cell is in the sample brain

[sample, centroid1, scale1] = tps_normalize(sample2);
% [cell, centroid2, scale2] = tps_normalize(cell);
% 
%[ref1, centroid3, scale3] = tps_normalize(ref);
% centroid4 = mean(regist);
 
[n,dim] = size(cell);
[m,dim] = size(ref);

reconcell = zeros(n,dim);
for i = 1:n
    for j = 1:dim
        K = tps_kernel(ref, cell(i,:));
        reconcell(i,j) = affine(1,j)+cell(i,:)*affine(2:end,j); 
        % apply the same transformation from the samlpe brain to the cell                                                      
        reconcell(i,j) = reconcell(i,j) + w(:,j)'*K;
%         reconcell(i,j) = reconcell(i,j) + w(:,j)'*K;
    end
end


reconcell = tps_denormalize(reconcell, centroid1, scale1); 

figure;
scatter3(ref(:,1),ref(:,2),ref(:,3),'b+');
hold on;
scatter3(reconcell(:,1),reconcell(:,2),reconcell(:,3),'ro');

figure;
scatter3(ref(:,1),ref(:,2),ref(:,3),'b+');
hold on;
scatter3(cell(:,1),cell(:,2),cell(:,3),'ro');
% reconstruct the cell onto the ref brain

% %%
% [n,dim] = size(centroid2);
% recent = zeros(n,dim);
% [centd, centroid5, scale5] = tps_normalize(centroid2);
% for i = 1:n
%     for j = 1:dim
%         K2 = tps_kernel(ref1, centd);
%         recent(i,j) = affine(1,j)+centd(i,:)*affine(2:end,j);
%         recent(i,j) = recent(i,j) + w(:,j)'*K2;
%     end
% end
% 
% figure;
% scatter3(recent(:,1), recent(:,2),recent(:,3),'m*');
% hold on;
% scatter3(ref1(:,1),ref1(:,2),ref1(:,3),'b+');

%%
% for i = 1:n
% 
%     K = tps_kernel(ref, cell(i,:));
%     fx(i) = affine(1,1) +cell(i,:)*affine(2:end,1);
%     fy(i) = affine(1,2) + cell(i,:)*affine(2:end,2);
%     fz(i) = affine(1,3) + cell(i,:)*affine(2:end,3);
% % 
% 
%     fx(i) = fx(i)+ w(:,1)'*K;
%     fy(i) = fy(i) + w(:,2)'*K;
%     fz(i) = fz(i) + w(:,3)'*K;
% end
    
% K = tps_kernel(ref, cell);
% fx = affine(1,1) +cell*affine(2:end,1);
% fy = affine(1,2) + cell*affine(2:end,2);
% fz = affine(1,3) + cell*affine(2:end,3);


% 

% reconcell = [fx; fy; fz]';
% 
reconcell = tps_denormalize(reconcell,centroid3 , scale1/scale2); 
% ref = tps_denormalize(ref, centroid3, scale3);
% cell = tps_denormalize(cell, centroid2, scale2);

figure;
scatter3(ref(:,1),ref(:,2),ref(:,3),'b+');
hold on;
scatter3(reconcell(:,1),reconcell(:,2),reconcell(:,3),'ro');
% scatter3(centroid3(:,1), centroid3(:,2), centroid3(:,3), 'm*');
% scatter3(centroid4(:,1), centroid4(:,2), centroid4(:,3), 'y*');
% figure; 
% scatter3(sample2(:,1),sample2(:,2),sample2(:,3),'b+');
% hold on;
% scatter3(cell(:,1),cell(:,2),cell(:,3),'ro');
% scatter3(centroid1(:,1), centroid1(:,2), centroid1(:,3), 'm*');
% scatter3(centroid2(:,1), centroid2(:,2), centroid2(:,3), 'y*');
% hold on;
% plot3(cell(:,1), cell(:,2),cell(:,3), 'g*');
% hold on;
% plot3(sample5(:,1), sample5(:,2),sample5(:,3), 'm+');

% 
% %%
% fx = a1 + x*affine(2:end,1) + w(:,j)*norm(ref(i,:)-cell(i,:))




% K = tps_kernel(sample, sample);
% 
% Pn = [ones(n,1), sample]; % pi = (1, xi, yi, zi)
% 
% Ls = [K Pn; Pn' zeros(4,4)];
% param = Ls\[ref; zeros(4,3)]; 
% 
% w = param(1:n,:);
% affine = param(n+1:end, :);
 

% U = tps_kernel(ref,cell);
% % bestlam = findlam(sample, ref);
% % U = U + bestlam*U;
% 
% P = [ones(m,1), ref];
% L = [U,P];
% regist = L * param;












        