function [Loading Score Eig]=pca(X, numPC)
% Principal Conponent Analysis
[n,p]=size(X);
C=X'*X/(n-1);
[V D]=eig(C); % V = Eigen Vector, D= Eigenvalues
d=diag(D);
[Eig,Index]=sort(d,'descend');

if nargin < 2
    numPC = p;
end

for i = 1:numPC
    Loading(:,i) = V(:,Index(i))/norm(V(:,Index(i)));
end

Score=X*Loading;    