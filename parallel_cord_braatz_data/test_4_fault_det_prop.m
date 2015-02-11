%Here I am trying to do fault detection for other faults by using only 
%first 160 samples from NO data
% Things changed : avg calc for NO & sd data, 

clear all
clc
% Fault_PCA
%% Phase 1: off-line modeling
% step 1: input data
load d00_te.mat; % sampled data under normal condition
s=[1:22,42:52];                    % variable index/3 of columns 
ss=[161:960];                      % how many time samples
X0=d00_te(:,s);
 
% step 2: data normalization
[Xref,Xmean,Xstd]=zscore(X0);           % data normalization
[N,p] = size(Xref);


% step 3: principal component analysis
[Pref,Tref,Eref]=pca(Xref);              % Loading score Eigen     
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
load d00_te.mat;load d01_te.mat;load d02_te.mat;load d03_te.mat;    
load d04_te.mat;load d05_te.mat;load d06_te.mat;load d07_te.mat;
load d08_te.mat;load d09_te.mat;load d10_te.mat;load d11_te.mat;  
load d12_te.mat;load d13_te.mat;load d14_te.mat;load d15_te.mat;  
load d16_te.mat;load d17_te.mat;load d18_te.mat;load d19_te.mat;
load d20_te.mat;load d21_te.mat;
load d00_05sd_te.mat;load d00_1sd_te.mat;load d00_2sd_te.mat;

    for i=1:9
        eval(['X','=','d0',num2str(i),'_te(:,s)',';']); %change to _te(ss,s) if you want remove NO data
        eval(['Xcrt','=','autoscale(X,Xmean,Xstd)',';']);
        eval(['Xcrt',num2str(i),'=','Xcrt',';']);
    end
    
    for i=10:21
        eval(['X','=','d',num2str(i),'_te(:,s)',';']); %change to _te(ss,s) if you want remove NO data
        eval(['Xcrt','=','autoscale(X,Xmean,Xstd)',';']);
        eval(['Xcrt',num2str(i),'=','Xcrt',';']);
    end
    
    X=d00_05sd_te(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling 
    Xcrt05sd = Xcrt;
    X=d00_1sd_te(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling 
    Xcrt1sd = Xcrt;
    X=d00_2sd_te(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling 
    Xcrt2sd = Xcrt;
    
    X=d00_te(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling
    
    n = size(Xcrt,1);                         % # of rows in Xcrt
    nn = size(Xcrt1,1);                       % # of rows in Xcrt1 
    
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
    FaultyPCA = plot(1:size(Xcrt1,2), Xcrt1(740:780,:),'LineStyle','--',...
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
    FaultyPCA = plot(1:size(Xcrt3,2), Xcrt3(740:780,:),'LineStyle','--',...
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
    Tcrt05sd = Xcrt05sd*Pref(:,1:a);
    Tcrt1sd = Xcrt1sd*Pref(:,1:a);
    Tcrt2sd = Xcrt2sd*Pref(:,1:a);
    % for loop to find Tcrt for all faults 
    for i=1:21
        eval(['Tcrt',num2str(i),'=','Xcrt',num2str(i),'*','Pref(:, 1:a)',';']);
    end
    
% %step 4: monitoring statistic for fault # 1
    T2=zeros(n,1);
    Q=zeros(n,1);
for i=1:n
    T2(i)=Tcrt1(i,:)*inv(diag(Eref(1:a)))*Tcrt1(i,:)';
    Q(i)=Xcrt1(i,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')*Xcrt1(i,:)';
end

% For T2 and Q of other faults 
% for i=1:n
%     T2(i)=Tcrt21(i,:)*inv(diag(Eref(1:a)))*Tcrt21(i,:)';
%     Q(i)=Xcrt21(i,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')*Xcrt21(i,:)';
% end

% %step 5: monitoring profiles
    N0 = 160;
    figure()
    subplot(1,2,1)
    plot(1:N0,T2(1:N0),'b-')
    hold on
    plot(N0+1:n,T2(N0+1:n),'r-')
    xlabel('Samples Fault#1')
    ylabel('T^2')
    plot(1:n,repmat(T2c,1,n),'k-')
    hold off
    subplot(1,2,2)
    plot(1:N0,Q(1:N0),'b-')
    hold on
    plot(N0+1:n,Q(N0+1:n),'r-')
    xlabel('Samples Fault#1')
    ylabel('Q')
    plot(1:n,repmat(Qc,1,n),'k-')
    hold off

%% Multidimensional Visualization 
    figure()
    Normref = plot(1:size(Tcrt,2), Tcrt(120:160,:),'LineStyle','--',...
        'Marker',markers(1,1),'color',color_code(1,:),'linewidth',1);hold on
    FaultyPCA = plot(1:size(Tcrt1,2), Tcrt1(740:760,:),'LineStyle','--',...
        'Marker',markers(1,4),'color',color_code(2,:),'linewidth',1);hold on
    FaultyPCAadd = plot(1:size(Tcrt1,2), Tcrt1(761:780,:),'LineStyle','--',...
        'Marker',markers(1,4),'color',color_code(2,:),'linewidth',1);hold on
    NSGroup = hggroup; 
    hSGroup = hggroup; 
    hCGroup = hggroup;
    set(Normref,'Parent',NSGroup)
    set(FaultyPCAadd,'Parent',hSGroup)
    set(FaultyPCA,'Parent',hCGroup)
    % Include these hggroups in the legend:
    set(get(get(NSGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    set(get(get(hSGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    set(get(get(hCGroup,'Annotation'),'LegendInformation'),...
        'IconDisplayStyle','on'); 
    legend('Scores from Normal Operating Data (Samples 120 to 160)',...
        'Scores from Faulty Operating Data (Samples 900 to 940)')
    grid on
    axis([1 10 -20 10]);
    set(gca,'XTick',1:1:10);
    xlabel('Principal Component');
    ylabel('Scores');
    hold off

%% min,max,avg,std dev, CI calculations

mini = [];mini05sd =[];mini1sd =[];mini2sd =[]; mini_modified3=[];
maxi = []; maxi05sd=[];maxi1sd =[];maxi2sd =[];maxi_modified3=[];
avg_data = [];avg_data05sd = [];avg_data1sd = [];avg_data2sd = [];
avg_datamodified3 = []; stdPCmodified3 = [];
std_PC = [];std_PC05sd = [];std_PC1sd = [];std_PC2sd = [];

CI = [];

% for loop to find naming max for all faults 
 
    for i=1:21
        eval(['mini',num2str(i),'=','[]'';']);
        eval(['avg_data',num2str(i),'=','[]',';']);
        eval(['maxi',num2str(i),'=','[]'';']);
        eval(['std_PC',num2str(i),'=','[]'';']);
    end

% for loop to calc for all faults 
    for i=1:21                          % # of fault data sets
        for j=1:size(Tcrt,2)            % # of columns in Tcrt
        eval(['mini',num2str(i),'=','[','mini',num2str(i),';',...
            'min(Tcrt',num2str(i),'(:,j))]',';']);
        eval(['maxi',num2str(i),'=','[','maxi',num2str(i),';',...
            'max(Tcrt',num2str(i),'(:,j))]',';']);
        eval(['avg_data',num2str(i),'(j,1)','=',...
            'mean(Tcrt',num2str(i),'(:,j))',';']);
        eval(['std_PC',num2str(i),'=','[','std_PC',num2str(i),';',...
            'std(Tcrt',num2str(i),'(:,j))]',';']);
        end
    end

for i=1:size(Tcrt,2)                   % # of columns in Tcrt
    mini = [mini; min(Tcrt(:,i))];
    maxi = [maxi; max(Tcrt(:,i))];
    avg_data(i,1) = mean(Tcrt(1:160,i)); %only first 160 samples used for average 
    std_PC = [std_PC; std(Tcrt(:,i))]; %only first 160 samples used for average 
end

for i=1:size(Tcrt05sd,2)
    mini05sd = [mini05sd; min(Tcrt05sd(:,i))];
    maxi05sd = [maxi05sd; max(Tcrt05sd(:,i))];
    avg_data05sd(i,1) =  mean(Tcrt05sd(1:160,i));%only first 160 samples used for average 
    std_PC05sd = [std_PC05sd; std(Tcrt05sd(:,i))];%only first 160 samples used for average 
end

for i=1:size(Tcrt1sd,2)
    mini1sd = [mini1sd; min(Tcrt1sd(:,i))];
    maxi1sd = [maxi1sd; max(Tcrt1sd(:,i))];
    avg_data1sd(i,1) =  mean(Tcrt1sd(1:160,i));%only first 160 samples used for average 
    std_PC1sd = [std_PC1sd; std(Tcrt1sd(:,i))];%only first 160 samples used for average 
end

for i=1:size(Tcrt2sd,2)
    mini2sd = [mini2sd; min(Tcrt2sd(:,i))];
    maxi2sd = [maxi2sd; max(Tcrt2sd(:,i))];
    avg_data2sd(i,1) =  mean(Tcrt2sd(1:160,i));%only first 160 samples used for average 
    std_PC2sd = [std_PC2sd; std(Tcrt2sd(:,i))];%only first 160 samples used for average 
end

% z_score = input('Specify the Z Score: ');
z_score = [1.96;2.58]; % 95,99 CI

for i = 1:length(z_score)
    CI_up(:,i) = avg_data + z_score(i)*(std_PC./sqrt(160));%size(Tcrt,1)) only first 160 samples used for average 
    CI_low(:,i) = avg_data - z_score(i)*(std_PC./sqrt(160));%only first 160 samples used for average 
end

figure()
plot(1:size(Tcrt,2),Tcrt(:,:))
figure()
plot(1:size(Tcrt,2),Tcrt1(:,:))
figure()
plot(1:size(Tcrt,2),Tcrt3(:,:))


%% Fault Detection Visualization 
figure()
subplot(2,1,1);
plot([1:size(Tcrt,2)], avg_data, 'bo-','linewidth',1.15); hold on
plot([1:size(Tcrt1,2)], avg_data1, 'r*-','linewidth',1.15);hold on
plot([1:size(Tcrt2,2)], avg_data2, 'g>-','linewidth',1.15);hold on
plot([1:size(Tcrt3,2)], avg_data3, 'md-','linewidth',1.15);hold on;

plot([1:size(Tcrt,2)], CI_up(:,2),'k-.','linewidth',1);%99 CI
plot([1:size(Tcrt,2)], CI_low(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_up(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_low(:,2),'k-.','linewidth',1);%99 CI
xlabel('Principal Components')
ylabel('Mean and CI')
legend ('Normal Operating Data','Fault # 1','Fault # 2',...
    'Fault # 3','99% CI','95% CI');
grid on

subplot(2,1,2);
plot([1:size(Tcrt,2)], avg_data, 'bo-','linewidth',1.15); hold on
plot([1:size(Tcrt1,2)], avg_data1, 'r*-','linewidth',1.15);hold on
plot([1:size(Tcrt2,2)], avg_data2, 'g>-','linewidth',1.15);hold on
plot([1:size(Tcrt3,2)], avg_data3, 'md-','linewidth',1.15);hold on;

plot([1:size(Tcrt,2)], CI_up(:,2),'k-.','linewidth',1);%99 CI
plot([1:size(Tcrt,2)], CI_low(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_up(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_low(:,2),'k-.','linewidth',1);%99 CI
grid on
axis([1 10 -0.7 1.2])
xlabel('Principal Components')
ylabel('Mean and CI')
hold off;

%% Plot Std dev
figure()
subplot(2,1,1);
plot([1:size(Tcrt,2)], std_PC, 'b-','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',8); hold on
plot([1:size(Tcrt1,2)], std_PC1, 'rv','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',8); hold on
plot([1:size(Tcrt05sd,2)], std_PC05sd, 'gp','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',8); hold on
plot([1:size(Tcrt1sd,2)], std_PC1sd, 'm-','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1],'MarkerSize',8); hold on
plot([1:size(Tcrt2sd,2)], std_PC2sd, 'g--','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1],'MarkerSize',8); hold on
plot([1:size(Tcrt3,2)], std_PC3, 'r^','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',8); hold on

xlabel('Principal Components')
ylabel('Standard Deviation')
legend ('Normal Operating Data','Fault # 1 Operating Data','Fault # 05sd Operating Data',...
    'Fault # 1sd Operating Data','Fault # 2sd Operating Data','Fault # 3 Operating Data');
grid on

subplot(2,1,2);
plot([1:size(Tcrt,2)], std_PC, 'b-','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',8); hold on
plot([1:size(Tcrt1,2)], std_PC1, 'rv','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',8); hold on
plot([1:size(Tcrt05sd,2)], std_PC05sd, 'gp','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',8); hold on
plot([1:size(Tcrt1sd,2)], std_PC1sd, 'm-','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1],'MarkerSize',8); hold on
plot([1:size(Tcrt2sd,2)], std_PC2sd, 'g--','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1],'MarkerSize',8); hold on
plot([1:size(Tcrt3,2)], std_PC3, 'r^','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',8); hold on

grid on
axis([1 10 1 3])
xlabel('Principal Components')
ylabel('Standard Deviation')
legend ('Normal Operating Data','Fault # 1 Operating Data','Fault # 05sd Operating Data',...
    'Fault # 1sd Operating Data','Fault # 2sd Operating Data','Fault # 3 Operating Data');

hold off;

%% Area/Centroid of Trapezoid calculation
%Center of area or center of mass for a uniform lamina
centroid = [];centroid05sd = [];centroid1sd = [];centroid2sd = [];
centroid_modified3 = [];area_modified3 = [];
area = []; area05sd = []; area1sd = [];area2sd = [];
    for i=1:21
        eval(['centroid',num2str(i),'=','[]',';']);
        eval(['area',num2str(i),'=','[]',';']);
    end

    for i=1:21
        for j=1:size(Tcrt,2)-1
        eval(['bb','=','(max(Tcrt',num2str(i),'(:,j))','-','min(Tcrt',num2str(i),'(:,j)))',';']);
        eval(['aa','=','(max(Tcrt',num2str(i),'(:,j+1))','-','min(Tcrt',num2str(i),'(:,j+1)))',';']);
        eval(['centroid',num2str(i),'(j)','=','(bb+2*aa)','/','(bb+aa)','/','3',';']);
        eval(['area',num2str(i),'(j)','=','(bb+aa)','/','2',';']);
        end
    end    
    
for i=1:size(Tcrt,2)-1
    bb = max(Tcrt(:,i)) - min(Tcrt(:,i));
    aa = max(Tcrt(:,i + 1)) - min(Tcrt(:,i + 1));
    centroid(i) = (bb + 2*aa)/(bb + aa)/3;
    area(i) = (bb + aa)/2;
end

for i=1:size(Tcrt05sd,2)-1
    bb = max(Tcrt05sd(:,i)) - min(Tcrt05sd(:,i));
    aa = max(Tcrt05sd(:,i + 1)) - min(Tcrt05sd(:,i + 1));
    centroid05sd(i) = (bb + 2*aa)/(bb + aa)/3;
    area05sd(i) = (bb + aa)/2;
end
for i=1:size(Tcrt1sd,2)-1
    bb = max(Tcrt1sd(:,i)) - min(Tcrt1sd(:,i));
    aa = max(Tcrt1sd(:,i + 1)) - min(Tcrt1sd(:,i + 1));
    centroid1sd(i) = (bb + 2*aa)/(bb + aa)/3;
    area1sd(i) = (bb + aa)/2;
end
for i=1:size(Tcrt2sd,2)-1
    bb = max(Tcrt2sd(:,i)) - min(Tcrt2sd(:,i));
    aa = max(Tcrt2sd(:,i + 1)) - min(Tcrt2sd(:,i + 1));
    centroid2sd(i) = (bb + 2*aa)/(bb + aa)/3;
    area2sd(i) = (bb + aa)/2;
end

%% Centroid plotting 
figure()
subplot(2,1,1);
aa = (1:size(Tcrt,2))-0.5;    
centroidaxis = aa(1,2:size(aa,2)); % used 2:size..to omit first x axis label
plot(centroidaxis,centroid,'b*');hold on

plot([1:size(Tcrt,2)], maxi, 'b-');hold on
plot([1:size(Tcrt1,2)], maxi1, 'r-');hold on
plot([1:size(Tcrt2,2)], maxi2, 'g-');hold on
plot([1:size(Tcrt3,2)], maxi3, 'm-');hold on

plot([1:size(Tcrt,2)], mini, 'b-'); hold on
plot([1:size(Tcrt1,2)], mini1, 'r-');hold on
plot([1:size(Tcrt2,2)], mini2, 'g-');hold on
plot([1:size(Tcrt3,2)], mini3, 'm-');hold on

plot(centroidaxis,centroid1,'r*');hold on
plot(centroidaxis,centroid2,'g*');hold on
plot(centroidaxis,centroid3,'m*');hold on
xlabel('Principal Components')
ylabel('Centroid & Min/Max of Principal components')
legend ('Centroid NO','Normal Operating data','Fault # 1','Fault # 2',...
    'Fault # 3');
grid on
%hold off;

subplot(2,1,2);
plot(centroidaxis,centroid,'b*');hold on
plot(centroidaxis,centroid1,'r*');hold on
plot(centroidaxis,centroid2,'g*');hold on
plot(centroidaxis,centroid3,'m*');hold on
axis([1 10 0.35 0.6])
grid on
hold off;

%% Trapezoid Area plotting 

figure()
subplot(2,1,1);
aa = (1:size(Tcrt,2))-0.5;    
centroidaxis = aa(1,2:size(aa,2)); % used 2:size..to omit first x axis label
plot(centroidaxis,area,'b*');hold on

plot([1:size(Tcrt,2)], maxi, 'b-');hold on
plot([1:size(Tcrt1,2)], maxi1, 'r-');hold on
plot([1:size(Tcrt2,2)], maxi2, 'g-');hold on
plot([1:size(Tcrt3,2)], maxi3, 'm-');hold on

plot([1:size(Tcrt,2)], mini, 'b-'); hold on
plot([1:size(Tcrt1,2)], mini1, 'r-');hold on
plot([1:size(Tcrt2,2)], mini2, 'g-');hold on
plot([1:size(Tcrt3,2)], mini3, 'm-');hold on

plot(centroidaxis,area1,'r*');hold on
plot(centroidaxis,area2,'g*');hold on
plot(centroidaxis,area3,'m*');hold on
xlabel('Principal Components')
ylabel('Area & Min/Max of Principal components')
legend ('Area NO','Normal Operating data','Fault # 1','Fault # 2',...
    'Fault # 3');
grid on

subplot(2,1,2);
plot(centroidaxis,area,'b*');hold on
plot(centroidaxis,area1,'r*');hold on
plot(centroidaxis,area2,'g*');hold on
plot(centroidaxis,area3,'m*');hold on
grid on
hold off;

%% Samples at different time
    
    figure()
    subplot(1,2,1);
    plot(1:size(Tcrt1,2), Tcrt1(80,:),'k-','linewidth',1.15);hold on
    plot(1:size(Tcrt1,2), Tcrt1(180,:),'b-','linewidth',1.15);hold on
    plot(1:size(Tcrt1,2), Tcrt1(280,:),'r-','linewidth',1.15);hold on
    plot(1:size(Tcrt1,2), Tcrt1(380,:),'g-','linewidth',1.15);hold on
    plot(1:size(Tcrt1,2), Tcrt1(480,:),'m-','linewidth',1.15);hold on
    plot(1:size(Tcrt1,2), Tcrt1(580,:),'c-','linewidth',1.15);hold on
    title('Samples 80,180,280,380,480 and 580');
    xlabel('Principal Components')
    ylabel('Score')
    legend ('Sample 80','Sample 180','Sample 280','Sample 380','Sample 480',...
        'Sample 580');
    axis([1 10 -25 20])
    grid on
    hold off;
    subplot(1,2,2);
    plot(1:size(Tcrt1,2), Tcrt1(480,:),'m-','linewidth',1.15);hold on
    plot(1:size(Tcrt1,2), Tcrt1(580,:),'c-','linewidth',1.15);hold on
    title('Samples 480 and 580');
    xlabel('Principal Components')
    ylabel('Score')
    legend ('Sample 480','Sample 580');
    axis([1 10 -20 10])
    grid on
    hold off;
    
%% Contribution Plots
SampleNo = 480;        % the sample No. you use to generate contribution plots

% contribution plots for Hotelling's T2
ConT1 = (Xcrt1(SampleNo, :)*Pref(:, 1:a)*sqrt(inv(diag(Eref(1:a))))*Pref(:, 1:a)').^2;
figure()
subplot(2,2,1)
bar(ConT1,'r');
xlabel('Process Variables'); ylabel('Contribution to T2');
legend ('Sample 480');
axis([0 34 -2 140])
grid on

% contribution plots for Q
ConQ1 = (Xcrt1(SampleNo,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')).^2;
subplot(2,2,3)
bar(ConQ1,'r');xlabel('Process Variables'); ylabel('Contribution to Q');
legend ('Sample 480');
axis([0 34 -2 25])
grid on

SampleNo = 580;     % the sample No. you use to generate contribution plots

% contribution plots for Hotelling's T2
ConT2 = (Xcrt1(SampleNo, :)*Pref(:, 1:a)*sqrt(inv(diag(Eref(1:a))))*Pref(:, 1:a)').^2;
subplot(2,2,2)
bar(ConT2,'b');xlabel('Process Variables'); ylabel('Contribution to T2');
legend ('Sample 580');
axis([0 34 -2 140])
grid on

% contribution plots for Q
ConQ2 = (Xcrt1(SampleNo,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')).^2;
subplot(2,2,4)
bar(ConQ2,'b');xlabel('Process Variables'); ylabel('Contribution to Q');
legend ('Sample 580');
axis([0 34 -2 25])
grid on

%% Plotting NO,F#1,F#3 and 0.5,1 & 2 std dev data with centroids/area

figure()
subplot(2,1,1);
plot([1:size(Tcrt,2)], avg_data, 'bo-','linewidth',1.15); hold on
plot([1:size(Tcrt1,2)], avg_data1, 'r*-','linewidth',1.15);hold on
plot([1:size(Tcrt3,2)], avg_data3, 'm*-','linewidth',1.15);hold on
plot([1:size(Tcrt05sd,2)], avg_data05sd, 'g>-','linewidth',1.15);hold on
plot([1:size(Tcrt1sd,2)], avg_data1sd, 'gd-','linewidth',1.15);hold on;
plot([1:size(Tcrt2sd,2)], avg_data2sd, 'cd-','linewidth',1.15);hold on;


plot([1:size(Tcrt,2)], CI_up(:,2),'k-.','linewidth',1);%99 CI
plot([1:size(Tcrt,2)], CI_low(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_up(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_low(:,2),'k-.','linewidth',1);%99 CI
xlabel('Principal Components')
ylabel('Mean and CI')
legend ('Normal Operating Data','Fault # 1','Fault # 3','0.5 sd',...
    '1 sd','2 sd','99% CI','95% CI');
grid on

subplot(2,1,2);
plot([1:size(Tcrt,2)], avg_data, 'bo-','linewidth',1.15); hold on
plot([1:size(Tcrt1,2)], avg_data1, 'r*-','linewidth',1.15);hold on
plot([1:size(Tcrt3,2)], avg_data3, 'm*-','linewidth',1.15);hold on
plot([1:size(Tcrt05sd,2)], avg_data05sd, 'g>-','linewidth',1.15);hold on
plot([1:size(Tcrt1sd,2)], avg_data1sd, 'gd-','linewidth',1.15);hold on;
plot([1:size(Tcrt2sd,2)], avg_data2sd, 'cd-','linewidth',1.15);hold on;

plot([1:size(Tcrt,2)], CI_up(:,2),'k-.','linewidth',1);%99 CI
plot([1:size(Tcrt,2)], CI_low(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_up(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_low(:,2),'k-.','linewidth',1);%99 CI
grid on
axis([1 10 -0.3 0.3])
xlabel('Principal Components')
ylabel('Mean and CI')
hold off;

% Centroid plotting 
figure()
subplot(2,1,1);
aa = (1:size(Tcrt,2))-0.5;    
centroidaxis = aa(1,2:size(aa,2)); % used 2:size..to omit first x axis label
plot(centroidaxis,centroid,'b*-');hold on
plot(centroidaxis,centroid1,'r*');hold on
plot(centroidaxis,centroid1,'m+');hold on
plot(centroidaxis,centroid05sd,'g*');hold on
plot(centroidaxis,centroid1sd,'m*');hold on
plot(centroidaxis,centroid2sd,'c*-');hold on

plot([1:size(Tcrt,2)], maxi, 'b-');hold on
plot([1:size(Tcrt1,2)], maxi1, 'r-');hold on
plot([1:size(Tcrt3,2)], maxi3, 'k-');hold on
plot([1:size(Tcrt05sd,2)], maxi05sd, 'g-');hold on
plot([1:size(Tcrt1sd,2)], maxi1sd, 'm-');hold on
plot([1:size(Tcrt2sd,2)], maxi2sd, 'c-');hold on

plot([1:size(Tcrt,2)], mini, 'b-'); hold on
plot([1:size(Tcrt1,2)], mini1, 'r-');hold on
plot([1:size(Tcrt3,2)], mini3, 'k-');hold on
plot([1:size(Tcrt05sd,2)], mini05sd, 'g-');hold on
plot([1:size(Tcrt1sd,2)], mini1sd, 'm-');hold on
plot([1:size(Tcrt2sd,2)], mini2sd, 'c-');hold on

xlabel('Principal Components')
ylabel('Centroid & Min/Max of Principal components')
legend ('Centroid NO','Centroid Fault # 1','Centroid Fault # 3',...
    'Centroid 0.5 sd','Centroid 1 sd','Centroid 2 sd');
grid on

subplot(2,1,2);
plot(centroidaxis,centroid,'b*-');hold on
plot(centroidaxis,centroid1,'r*');hold on
plot(centroidaxis,centroid3,'m+');hold on
plot(centroidaxis,centroid05sd,'g*');hold on
plot(centroidaxis,centroid1sd,'m*');hold on
plot(centroidaxis,centroid2sd,'c*-');hold on
axis([1 10 0.35 0.6])
grid on
hold off;

% Trapezoid Area plotting 

figure()
subplot(2,1,1);
aa = (1:size(Tcrt,2))-0.5;    
centroidaxis = aa(1,2:size(aa,2)); % used 2:size..to omit first x axis label

plot(centroidaxis,area,'b*');hold on
plot(centroidaxis,area1,'r*');hold on
plot(centroidaxis,area3,'m+');hold on
plot(centroidaxis,area05sd,'g*');hold on
plot(centroidaxis,area1sd,'m*');hold on
plot(centroidaxis,area2sd,'c*');hold on

plot([1:size(Tcrt,2)], maxi, 'b-');hold on
plot([1:size(Tcrt1,2)], maxi1, 'r-');hold on
plot([1:size(Tcrt3,2)], maxi3, 'k-');hold on
plot([1:size(Tcrt05sd,2)], maxi05sd, 'g-');hold on
plot([1:size(Tcrt1sd,2)], maxi1sd, 'm-');hold on
plot([1:size(Tcrt2sd,2)], maxi2sd, 'c-');hold on

plot([1:size(Tcrt,2)], mini, 'b-'); hold on
plot([1:size(Tcrt1,2)], mini1, 'r-');hold on
plot([1:size(Tcrt3,2)], mini3, 'k-');hold on
plot([1:size(Tcrt05sd,2)], mini05sd, 'g-');hold on
plot([1:size(Tcrt1sd,2)], mini1sd, 'm-');hold on
plot([1:size(Tcrt2sd,2)], mini2sd, 'c-');hold on

xlabel('Principal Components')
ylabel('Area & Min/Max of Principal components')
legend ('area NO','area Fault # 1','area Fault # 3','area 0.5 sd'...
    ,'area 1 sd','area 2 sd');
grid on

subplot(2,1,2);
plot(centroidaxis,area,'b*');hold on
plot(centroidaxis,area1,'r*');hold on
plot(centroidaxis,area3,'m+');hold on
plot(centroidaxis,area05sd,'g*');hold on
plot(centroidaxis,area1sd,'m*');hold on
plot(centroidaxis,area2sd,'c*');hold on
axis([1 10 5 15])
grid on
hold off;

%% Plot centroids for a small sample set

centroidsample = [];
centroidsample1 = [];
centroidsample2 = [];
centroidsample3 = [];


for i=1:size(Tcrt,2)-1
    j = 560:570;                % sample numbers should be more than one
                                %If you change j here, change below also
    bb = max(Tcrt(j,i)) - min(Tcrt(j,i));
    aa = max(Tcrt(j,i + 1)) - min(Tcrt(j,i + 1));
    centroidsample(i) = (bb + 2*aa)/(bb + aa)/3;

    bb = max(Tcrt1(j,i)) - min(Tcrt1(j,i));
    aa = max(Tcrt1(j,i + 1)) - min(Tcrt1(j,i + 1));
    centroidsample1(i) = (bb + 2*aa)/(bb + aa)/3;
    
    bb = max(Tcrt2(j,i)) - min(Tcrt2(j,i));
    aa = max(Tcrt2(j,i + 1)) - min(Tcrt2(j,i + 1));
    centroidsample2(i) = (bb + 2*aa)/(bb + aa)/3;
    
    bb = max(Tcrt3(j,i)) - min(Tcrt3(j,i));
    aa = max(Tcrt3(j,i + 1)) - min(Tcrt3(j,i + 1));
    centroidsample3(i) = (bb + 2*aa)/(bb + aa)/3;
end

minisample = [];minisample1 = [];minisample2 = [];minisample3 = [];
maxisample = [];maxisample1 = [];maxisample2 = [];maxisample3 = [];

for i=1:size(Tcrt,2)
    j=560:564;              %If j is changed here, change above too
    minisample = [minisample; min(Tcrt(j,i))];
    maxisample = [maxisample; max(Tcrt(j,i))];
  
    minisample1 = [minisample1; min(Tcrt1(j,i))];
    maxisample1 = [maxisample1; max(Tcrt1(j,i))];

    minisample2 = [minisample2; min(Tcrt2(j,i))];
    maxisample2 = [maxisample2; max(Tcrt2(j,i))];

    minisample3 = [minisample3; min(Tcrt3(j,i))];
    maxisample3 = [maxisample3; max(Tcrt3(j,i))];

end

figure()
subplot(3,1,1);
aa = (1:size(Tcrt,2))-0.5;    
centroidaxis = aa(1,2:size(aa,2)); % used 2:size..to omit first x axis label
plot(centroidaxis,centroidsample,'b*');hold on

plot([1:size(Tcrt,2)], maxisample, 'b-');hold on
plot([1:size(Tcrt1,2)], maxisample1, 'r-');hold on
plot([1:size(Tcrt2,2)], maxisample2, 'g-');hold on
plot([1:size(Tcrt3,2)], maxisample3, 'm-');hold on

plot([1:size(Tcrt,2)], minisample, 'b-'); hold on
plot([1:size(Tcrt1,2)], minisample1, 'r-');hold on
plot([1:size(Tcrt2,2)], minisample2, 'g-');hold on
plot([1:size(Tcrt3,2)], minisample3, 'm-');hold on

plot(centroidaxis,centroidsample1,'r*');hold on
plot(centroidaxis,centroidsample2,'g*');hold on
plot(centroidaxis,centroidsample3,'m*');hold on
xlabel('Principal Components')
ylabel('Centroid & Min/Max of Principal components')
legend ('Centroid NO','Normal Operating data','Fault # 1','Fault # 2',...
    'Fault # 3');
grid on

subplot(3,1,2);
plot(centroidaxis,centroidsample,'b*');hold on
plot(centroidaxis,centroidsample1,'r*');hold on
plot(centroidaxis,centroidsample2,'g*');hold on
plot(centroidaxis,centroidsample3,'m*');hold on
axis([1 10 0.35 0.6])
grid on

subplot(3,1,3);
aa = (1:size(Tcrt,2))-0.5;    
centroidaxis = aa(1,2:size(aa,2));
plot([1:size(Tcrt,2)], maxisample, 'b-');hold on
%plot([1:size(Tcrt1,2)], maxisample1, 'r-');hold on
plot([1:size(Tcrt,2)], minisample, 'b-'); hold on
%plot([1:size(Tcrt1,2)], minisample1, 'r-');hold on
plot(centroidaxis,centroidsample,'b*');hold on
%plot(centroidaxis,centroidsample1,'r*');hold on

grid on
hold off;

%% For Friday 7 Nov meeting Fault 3,9,15,21 Detection Visualization 
figure()
subplot(2,1,1);
plot([1:size(Tcrt,2)], avg_data, 'LineStyle','-',...
   'Marker',markers(1,1),'color',color_code(1,:),'linewidth',1.15); hold on
plot([1:size(Tcrt1,2)], avg_data1,'LineStyle','--',...
   'Marker',markers(1,2),'color',color_code(2,:),'linewidth',1.15);hold on
plot([1:size(Tcrt3,2)], avg_data3,'LineStyle','--',...
   'Marker',markers(1,4),'color',color_code(5,:),'linewidth',1.15);hold on
plot([1:size(Tcrt9,2)], avg_data9,'LineStyle','--',...
   'Marker',markers(1,10),'color',color_code(11,:),'linewidth',1.15);hold on
plot([1:size(Tcrt15,2)], avg_data15,'LineStyle','--',...
   'Marker',markers(1,16),'color',color_code(17,:),'linewidth',1.15);hold on
plot([1:size(Tcrt21,2)], avg_data21,'LineStyle','--',...
   'Marker',markers(1,22),'color',color_code(23,:),'linewidth',1.15);hold on

plot([1:size(Tcrt,2)], CI_up(:,2),'k-.','linewidth',1);%99 CI
plot([1:size(Tcrt,2)], CI_low(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_up(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_low(:,2),'k-.','linewidth',1);%99 CI
xlabel('Principal Components')
ylabel('Mean and CI')
legend ('Normal Operating Data','Fault # 1','Fault # 3',...
    'Fault # 9','Fault # 15','Fault # 21','99% CI','95% CI');
grid on

subplot(2,1,2);
plot([1:size(Tcrt,2)], avg_data, 'LineStyle','-',...
   'Marker',markers(1,1),'color',color_code(1,:),'linewidth',1.15); hold on
plot([1:size(Tcrt1,2)], avg_data1,'LineStyle','--',...
   'Marker',markers(1,2),'color',color_code(2,:),'linewidth',1.15);hold on
plot([1:size(Tcrt3,2)], avg_data3,'LineStyle','--',...
   'Marker',markers(1,4),'color',color_code(5,:),'linewidth',1.15);hold on
plot([1:size(Tcrt9,2)], avg_data9,'LineStyle','--',...
   'Marker',markers(1,10),'color',color_code(11,:),'linewidth',1.15);hold on
plot([1:size(Tcrt15,2)], avg_data15,'LineStyle','--',...
   'Marker',markers(1,16),'color',color_code(17,:),'linewidth',1.15);hold on
plot([1:size(Tcrt21,2)], avg_data21,'LineStyle','--',...
   'Marker',markers(1,22),'color',color_code(23,:),'linewidth',1.15);hold on

plot([1:size(Tcrt,2)], CI_up(:,2),'k-.','linewidth',1);%99 CI
plot([1:size(Tcrt,2)], CI_low(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_up(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_low(:,2),'k-.','linewidth',1);%99 CI
grid on
axis([1 10 -0.3 0.3])
xlabel('Principal Components')
ylabel('Mean and CI')
hold off;







