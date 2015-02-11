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

% step 3: principal component analysis
[Pref,Tref,Eref]=pca(Xref);             
getpercent = 0.7;
a = cpv(Eref,getpercent);                     % number of PCs 

% step4: thresholds 
alpha=0.99;                                        % confidence level
T2c=a*(N-1)/(N-a)*finv(alpha,a,N-a);  % control limit of T

% for i=1:3
%     c(i)=sum(Eref(a+1:p).^i);
% end
% h0=1-2*c(1)*c(3)/3/c(2)^2;
% ca=norminv(alpha,0,1);
% Qc=c(1)*(ca*sqrt(2*c(2)*(h0^2))/c(1)+1+c(2)*h0*(h0-1)/c(1)^2)^(1/h0);

Xe = Tref(:,1:a)*Pref(:,1:a)';
Eror = Xref - Xe;
for i = 1:N
    Q(i,1) = Eror(i,:)*Eror(i,:)';
end
m = mean(Q); v = var(Q);
g = v/m/2; h = 2*m^2/v;
Qc = g * chi2inv(alpha,h);                    % control limit of SPE

%% Phase 2:  on-line monitoring
% step 1: on-line sampled data
load d05_te.mat;    
X=d05_te(:,s);

% step 2: data normalization
Xcrt=autoscale(X,Xmean,Xstd);              % data scaling
n = size(Xcrt,1);

% step 3: pca transformation
Tcrt = Xcrt*Pref(:, 1:a);

% step 4: monitoring statistics
T2=zeros(n,1);
Q=zeros(n,1);
for i=1:n
    T2(i)=Tcrt(i,:)*inv(diag(Eref(1:a)))*Tcrt(i,:)';
    Q(i)=Xcrt(i,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')*Xcrt(i,:)';
end

% step 5: monitoring profiles
N0 = 160;
figure(1)
subplot(2,1,1)
plot(1:N0,T2(1:N0),'b-')
hold on
plot(N0+1:n,T2(N0+1:n),'r-')
xlabel('Samples')
ylabel('T^2')
plot(1:n,repmat(T2c,1,n),'k-')
hold off

subplot(2,1,2)
plot(1:N0,Q(1:N0),'b-')
hold on
plot(N0+1:n,Q(N0+1:n),'r-')
xlabel('Samples')
ylabel('Q')
plot(1:n,repmat(Qc,1,n),'k-')
hold off

%% Contribution Plots
SampleNo = 170;        % the sample No. you use to generate contribution plots
% contribution plots for Hotelling's T2
ConT = (Xcrt(SampleNo, :)*Pref(:, 1:a)*sqrt(inv(diag(Eref(1:a))))*Pref(:, 1:a)').^2;
figure(2)
subplot(2,1,1)
bar(ConT)

% contribution plots for Q
ConQ = (Xcrt(SampleNo,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')).^2;
subplot(2,1,2)
bar(ConQ)

%% Multidimensional Visualization 
figure(3)
hold on
for i = 1:N
    plot(1:a, Tref(i, 1:a), 'k');
end
grid on

n0 = 160;
for i = 1:n
    pause(0.2);
    if i > n0
        plot(1:a, Tcrt(i,:), 'r');
    else
        plot(1:a, Tcrt(i,:), 'b');
    end
end

hold off