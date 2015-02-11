clear
clc
% Fault_PCA
%% Phase 1: off-line modeling
% step 1: input data
load d00_te.mat;                                 % sampled data under normal condition
s=[1:22,42:52];                                  % variable index 
X0=d00_te(:,s);

% step 2: data normalization
[Xref,Xmean,Xstd]=zscore(X0);           % data normalization
[N,p] = size(Xref);

% plot(d00_te(:,), d00_te(4,20:30), 'r');
% hold on
% plot(20:30, X0(4,20:30), 'r');
% hold on
% plot(20:30, X0(804,20:30), 'b');
% hold on

% figure(7) % Plot TE process data
% normdata = 160;
% for i = 1:n
%     pause(0.2);
%     if i > normdata
%         plot(1:p, X0(i,:), 'r');
%     hold on
%     else
%         plot(1:p, X0(i,:), 'b');
%     hold on
%     end
% end
% plot(1:p, X0(500:504,:),'g');
% hold on
% hold off


% step 3: principal component analysis
[Pref,Tref,Eref]=pca(Xref);              % Loading score Eigen     
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

%% Multidimensional Visualization 
figure(2)
hold on
for i = 1:(N-956)
    plot(1:a, Tref(i, 1:a), 'k');
end
grid on

n0 = 160;
for i = 1:(n-956)
    pause(0.2);
    if i > n0
        plot(1:a, Tcrt(i,:), 'r');
    else
        plot(1:a, Tcrt(i,:), 'b');
    end
end
scatter3(Tcrt(:,1),Tcrt(:,2),Tcrt(:,3))
title('3 PC in x,y and z')
biplot(Tcrt(:,1:2))
plot(Tcrt(:,1),Tcrt(:,2));
hold on

plot3(Tcrt(:,1),Tcrt(:,2),Eror(:,1)); %plot Error
title('2 PC in x,y and error on z')

hold off