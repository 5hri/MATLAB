clear
clc

%% Phase 1: off-line modeling
% step 1: input data
load d00_te.mat;                                 % sampled data under normal condition
s=[1:22,42:52];                                  % variable index
X0=d00_te(:,s);

% step 2: data normalization
[Xref,Xmean,Xstd]=zscore(X0);           % data normalization
[N,p] = size(Xref);

figure
hold on
for i = 1:N
    plot(1:p, Xref(i, 1:p), 'k')
end
grid on

% step 3: Univariate 2-sigma visualization

XstdXref = std(Xref(:,1))
plot(1:N,Xref(:,1), 'k')
hold on
plot(1:N,repmat(2*XstdXref,1,N), 'r')
hold on
plot(1:N,repmat(-2*XstdXref,1,N), 'r')
hold off

load d05_te.mat;    %load faulty data
X=d05_te(:,s);

% step 2: data normalization
Xcrt=autoscale(X,Xmean,Xstd);              % data scaling
n = size(Xcrt,1);

n0 = 160;
for i = 1:n
    pause(0.1)
    if i > n0
        plot(1:p, Xcrt(i, :), 'r');
    else
        plot(1:p, Xcrt(i, :), 'b');
    end
end

hold off

% step 3: Univariate 2-sigma visualization

XstdXcrt = std(Xcrt(:,1))
plot(1:N,Xcrt(:,1), 'k')
hold on
plot(1:N,repmat(2*XstdXcrt,1,N), 'r')
hold on
plot(1:N,repmat(-2*XstdXcrt,1,N), 'r')
hold off