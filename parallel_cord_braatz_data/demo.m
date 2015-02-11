clear
clc

% Step 1: input data for analysis
load('depent_norm_3_5hr.mat');
 X0 = depent_norm_3_5hr(2:end, :);
 
 % Step 2: data normalization
 [X, Xmean, Xstd] = zscore(X0);
 [n,m] = size(X);
 
 % Step 3: principal component analysis
 [P, T, eigval] = pca(X);
 getpercent = 0.90;
a = cpv(eigval,getpercent);                                       % number of PCs
 
% Step 4: multidimensional visualization
 figure(1)
 hold on;
axis([0.6 a+0.4 -6 6]);
for i = 1:a
    plot([i,i], [-6,6], 'k')
end

 for i = 1:n
     plot(1:a, T(i, 1:a), 'b');
     pause(0.2);
 end
 hold off;