clear all
clc
% Fault_PCA
%% Phase 1: off-line modeling
% step 1: input data
load no_data.mat; % sampled data under normal condition
s=[1:30,32:33];                    % # of columns 
ss=size(no_data,1); %[161:960];                      % how many time samples
X0=no_data(:,s);
 
% step 2: data normalization
[Xref,Xmean,Xstd]=zscore(X0);           % data normalization, Xref is Zscore of X0
[N,p] = size(Xref);

% step 3: principal component analysis
[Pref,Tref,Eref]=pca(Xref);              % Loading score Eigenvalues     
getpercent = 0.7;
a = cpv(Eref,getpercent);                     % number of PCs 

% step4: thresholds 
alpha=0.99;                           % confidence level
T2c=a*(N-1)/(N-a)*finv(alpha,a,N-a);  % control limit of T

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

load f1_data.mat;load no_data.mat;load f3_data.mat;

    X=f1_data(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling 
    Xcrtf1 = Xcrt;
    X=f3_data(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling 
    Xcrtf3 = Xcrt;
    
    X=no_data(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling 
    Xcrtno = Xcrt;
    
    n = size(Xcrt,1);                         % # of rows in Xcrt
    nn = size(Xcrtf1,1);                       % # of rows in Xcrtf1 
    
%% Distinguishable color
    color_code=distinguishable_colors(30);
    markers=['o','x','+','*','v','d','s','^','<','>','p','h','.',...
                '+','*','o','x','^','<','h','.','>','p','v','d','s',...
                'o','x','+','*','v','d','s','^','<','>','p','h','.'];
    linestyles = cellstr(char('-','--','-.',':','-',':','-.','--','-',...
                ':','-',':','-.','--','-',':','-.','--','-',':','-.'));
%     d = reshape(color_code,[1 size(color_code)]);
%     figure
%     image(d);

%% Parallel Coord visualtion of faulty data
    figure()
    subplot(2,1,1)
    NormPCA = plot(1:size(Xcrt,2), Xcrt(120:160,:),'LineStyle','-',...
        'Marker',markers(1,1),'color',color_code(1,:));hold on
    FaultyPCA = plot(1:size(Xcrtf1,2), Xcrtf1(740:780,:),'LineStyle','--',...
        'Marker',markers(1,4),'color',color_code(2,:));
    hSGroup = hggroup;
    hCGroup = hggroup;
    set(NormPCA,'Parent',hSGroup)
    set(FaultyPCA,'Parent',hCGroup)
    % Include these hggroups in the legend:
    set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    legend('Normal Operating Data (Samples 120 to 160)',...
            'Faulty Operating Data #1 (Samples 900 to 940)')
        %740:780 beacuse 900-160=740, 940-160=780,160 is the NO samples in
        %faulty data 
    grid on
    axis([1 33 -12 20]);
    set(gca,'XTick',1:1:33);
    xlabel('Process Variables (Fault #1)');
    ylabel('Standardized operating data at different time step');
    
    subplot(2,1,2)
    NormPCA = plot(1:size(Xcrt,2), Xcrt(120:160,:),'LineStyle','-',...
                    'Marker',markers(1,1),'color',color_code(1,:));hold on
    FaultyPCA = plot(1:size(Xcrtf3,2), Xcrtf3(740:780,:),'LineStyle','--',...
                    'Marker',markers(1,4),'color',color_code(5,:));
    hSGroup = hggroup;
    hCGroup = hggroup;
    set(NormPCA,'Parent',hSGroup)
    set(FaultyPCA,'Parent',hCGroup)
    % Include these hggroups in the legend:
    set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    legend('Normal Operating Data (Samples 120 to 160)',...
            'Faulty Operating Data #3 (Samples 900 to 940)')
    grid on
    axis([1 33 -4 4]);
    set(gca,'XTick',1:1:33);
    xlabel('Process Variables (Fault #3)');
    ylabel('Standardized operating data at different time step');

%% PCA Related calculations

% % step 3: pca transformation
    Tcrt  = Xcrt*Pref(:, 1:a);
   
    Tcrtf1 = Xcrtf1*Pref(:,1:a);
    Tcrtf3 = Xcrtf3*Pref(:,1:a);
    Tcrtno = Xcrtno*Pref(:,1:a);
    
 
% %step 4: monitoring statistic for fault # 1
    T2=zeros(n,1);
    Q=zeros(n,1);
for i=1:n
    T2(i)=Tcrtf1(i,:)*inv(diag(Eref(1:a)))*Tcrtf1(i,:)';
    Q(i)=Xcrtf1(i,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')*Xcrtf1(i,:)';
end

% % % For T2 and Q of other faults 
% for i=1:n
%     T2(i)=Tcrtf3(i,:)*inv(diag(Eref(1:a)))*Tcrtf3(i,:)';
%     Q(i)=Xcrtf3(i,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')*Xcrtf3(i,:)';
% end

% %step 5: monitoring profiles
    N0 = 160;
    figure()
    subplot(1,2,1)
    plot(1:N0,T2(1:N0),'b-')
    hold on
    plot(N0+1:n,T2(N0+1:n),'ro','MarkerSize',2)
    xlabel('Samples Fault#1')
    ylabel('T^2')
    plot(1:n,repmat(T2c,1,n),'k-')
    hold off
    subplot(1,2,2)
    plot(1:N0,Q(1:N0),'b-')
    hold on
    plot(N0+1:n,Q(N0+1:n),'ro','MarkerSize',2)
    xlabel('Samples Fault#1')
    ylabel('Q')
    plot(1:n,repmat(Qc,1,n),'k-')
    hold off
    
% % For T2 and Q of other fault # 3
for i=1:n
    T2(i)=Tcrtf3(i,:)*inv(diag(Eref(1:a)))*Tcrtf3(i,:)';
    Q(i)=Xcrtf3(i,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')*Xcrtf3(i,:)';
end

% %step 5: monitoring profiles
    N0 = 160;
    figure()
    subplot(1,2,1)
    plot(1:N0,T2(1:N0),'b-')
    hold on
    plot(N0+1:n,T2(N0+1:n),'ro','MarkerSize',2)
    xlabel('Samples Fault#1')
    ylabel('T^2')
    plot(1:n,repmat(T2c,1,n),'k-')
    hold off
    subplot(1,2,2)
    plot(1:N0,Q(1:N0),'b-')
    hold on
    plot(N0+1:n,Q(N0+1:n),'ro','MarkerSize',2)
    xlabel('Samples Fault#1')
    ylabel('Q')
    plot(1:n,repmat(Qc,1,n),'k-')
    hold off

%% plot most varying variables in normalized data in fault # 3
temp = 2100:2500; % specify the samples to plot 


figure;
for i=1:32
plot(1:size(temp,2),repmat(mean(Xcrtno(temp,i))+1.96*std(Xcrtno(temp,i)),1,size(temp,2)),'k');hold on
plot(1:size(temp,2),repmat(mean(Xcrtno(temp,i))-1.96*std(Xcrtno(temp,i)),1,size(temp,2)),'k');hold on
plot(Xcrtno(temp,i),'b.'); hold on
plot(Xcrtf3(temp,i),'ro'); hold on
xlabel(i)
hold off
pause(2);

end

subplot(5,2,1)    
plot(1:size(temp,2),repmat(mean(Xcrtno(temp,7))+1.96*std(Xcrtno(temp,7)),1,size(temp,2)),'k');hold on
plot(1:size(temp,2),repmat(mean(Xcrtno(temp,7))-1.96*std(Xcrtno(temp,7)),1,size(temp,2)),'k');hold on
plot(Xcrtno(temp,7),'bo'); hold on
plot(Xcrtf3(temp,7),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('Reactor Pressure')

subplot(5,2,2)

plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,11))+1.96*std(Xcrtno(temp,11)),1,size(Xcrtno,1)),'k');hold on
plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,11))-1.96*std(Xcrtno(temp,11)),1,size(Xcrtno,1)),'k');hold on
plot(Xcrtno(temp,11),'bo'); hold on
plot(Xcrtf3(temp,11),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('Product Separator Temperature')

subplot(5,2,3)

plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,13))+1.96*std(Xcrtno(temp,13)),1,size(Xcrtno,1)),'k');hold on
plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,13))-1.96*std(Xcrtno(temp,13)),1,size(Xcrtno,1)),'k');hold on
plot(Xcrtno(temp,13),'bo'); hold on
plot(Xcrtf3(temp,13),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('Product Separator Pressure')

subplot(5,2,4)

plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,16))+1.96*std(Xcrtno(temp,16)),1,size(Xcrtno,1)),'k');hold on
plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,16))-1.96*std(Xcrtno(temp,16)),1,size(Xcrtno,1)),'k');hold on
plot(Xcrtno(temp,16),'bo'); hold on
plot(Xcrtf3(temp,16),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('Stripper Pressure')

subplot(5,2,5)

plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,18))+1.96*std(Xcrtno(temp,18)),1,size(Xcrtno,1)),'k');hold on
plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,18))-1.96*std(Xcrtno(temp,18)),1,size(Xcrtno,1)),'k');hold on
plot(Xcrtno(temp,18),'bo'); hold on
plot(Xcrtf3(temp,18),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('Stripper Temperature')

subplot(5,2,6)

plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,19))+1.96*std(Xcrtno(temp,19)),1,size(Xcrtno,1)),'k');hold on
plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,19))-1.96*std(Xcrtno(temp,19)),1,size(Xcrtno,1)),'k');hold on
plot(Xcrtno(temp,19),'bo'); hold on
plot(Xcrtf3(temp,19),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('Stripper steam flow')

subplot(5,2,7)

plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,20))+1.96*std(Xcrtno(temp,20)),1,size(Xcrtno,1)),'k');hold on
plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,20))-1.96*std(Xcrtno(temp,20)),1,size(Xcrtno,1)),'k');hold on
plot(Xcrtno(temp,20),'bo'); hold on
plot(Xcrtf3(temp,20),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('Compressor Work')


subplot(5,2,8)

plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,24))+1.96*std(Xcrtno(temp,24)),1,size(Xcrtno,1)),'k');hold on
plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,24))-1.96*std(Xcrtno(temp,24)),1,size(Xcrtno,1)),'k');hold on
plot(Xcrtno(temp,24),'bo'); hold on
plot(Xcrtf3(temp,24),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('E feed flow valve')


subplot(5,2,9)

plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,27))+1.96*std(Xcrtno(temp,27)),1,size(Xcrtno,1)),'k');hold on
plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,27))-1.96*std(Xcrtno(temp,27)),1,size(Xcrtno,1)),'k');hold on
plot(Xcrtno(temp,27),'bo'); hold on
plot(Xcrtf3(temp,27),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('Compressor recycle valve')


subplot(5,2,10)

plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,31))+1.96*std(Xcrtno(temp,31)),1,size(Xcrtno,1)),'k');hold on
plot(1:size(Xcrtno,1),repmat(mean(Xcrtno(temp,31))-1.96*std(Xcrtno(temp,31)),1,size(Xcrtno,1)),'k');hold on
plot(Xcrtno(temp,31),'bo'); hold on
plot(Xcrtf3(temp,31),'ro'); hold on

axis([1 size(temp,2) -3 3]);
xlabel('Stripper Steam valve')
