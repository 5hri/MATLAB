function Y=autoscale(X,Xmean,Xstd)
% Scale the data X using mean Xmean, and standard deviation Xstd
n=size(X,1);
Y=(X-repmat(Xmean,n,1))./repmat(Xstd,n,1);
