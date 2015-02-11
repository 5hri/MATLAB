clear all
clc
% Fault_PCA
%% Phase 1: off-line modeling
% step 1: input data
load d00_te.mat; % sampled daata under normal condition
s=[1:22,42:52];                    % variable index 
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
load d000_te.mat;load d005_te.mat;

    X=d01_te(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling
    Xcrt1 = Xcrt;
    X=d02_te(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling
    Xcrt2 = Xcrt;
    X=d03_te(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling
    Xcrt3 = Xcrt;
    X=d005_te(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling 
    Xcrt005 = Xcrt;
    X=d00_te(:,s);
    Xcrt=autoscale(X,Xmean,Xstd);             % data scaling
    
    n = size(Xcrt,1);
    
%% Parallel Coord visualtion of faulty data
    figure
    NormPCA = plot(1:size(Xcrt1,2), Xcrt1(120:160,:),'b');hold on
    FaultyPCA = plot(1:size(Xcrt1,2), Xcrt1(900:940,:),'r--');
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
            'Faulty Operating Data (Samples 900 to 940)')
    grid on
    axis([1 33 -10 20]);
    set(gca,'XTick',1:1:33);
    xlabel('Process Variables');
    ylabel('Standardized operating data at different time step');
    
% step 3: pca transformation
    Tcrt  = Xcrt*Pref(:, 1:a);
    Tcrt1 = Xcrt1*Pref(:, 1:a);
    Tcrt2 = Xcrt2*Pref(:, 1:a);
    Tcrt3 = Xcrt3*Pref(:, 1:a);
    Tcrt005 = Xcrt005*Pref(:,1:a);
    
% step 4: monitoring statistic
    T2=zeros(n,1);
    Q=zeros(n,1);
for i=1:n
    T2(i)=Tcrt1(i,:)*inv(diag(Eref(1:a)))*Tcrt1(i,:)';
    Q(i)=Xcrt1(i,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')*Xcrt1(i,:)';
end

% step 5: monitoring profiles
    N0 = 160;
    figure()
    subplot(1,2,1)
    plot(1:N0,T2(1:N0),'b-')
    hold on
    plot(N0+1:n,T2(N0+1:n),'r-')
    xlabel('Samples')
    ylabel('T^2')
    plot(1:n,repmat(T2c,1,n),'k-')
    hold off
    subplot(1,2,2)
    plot(1:N0,Q(1:N0),'b-')
    hold on
    plot(N0+1:n,Q(N0+1:n),'r-')
    xlabel('Samples')
    ylabel('Q')
    plot(1:n,repmat(Qc,1,n),'k-')
    hold off

%% Multidimensional Visualization 
    
    Normref = plot(1:size(Tcrt1,2), Tcrt1(120:160,:),'b','linewidth',1);hold on
    FaultyPCA = plot(1:size(Tcrt1,2), Tcrt1(900:920,:),'r--','linewidth',1);hold on
    FaultyPCAadd = plot(1:size(Tcrt1,2), Tcrt1(921:940,:),'r--','linewidth',1);hold on
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

%% Fault Detection Visualization 

mini = [];mini1 = [];mini2 = [];mini3 = [];
maxi = [];maxi1 = [];maxi2 = [];maxi3 = [];mini005 = [];maxi005=[];
avg_data = [];avg_data1 = [];avg_data2 = [];avg_data3 = [];avg_data005 = [];
std_PC = [];std_PC1 = [];std_PC2 = [];std_PC3 = [];std_PC005 = [];CI = [];

for i=1:size(Tcrt,2)
    mini = [mini; min(Tcrt(:,i))];
    maxi = [maxi; max(Tcrt(:,i))];
    avg_data(i,1) = 0.5*(min(Tcrt(:,i)) + max(Tcrt(:,i)));
    std_PC = [std_PC; std(Tcrt(:,i))];
   
end

for i=1:size(Tcrt1,2)
    mini1 = [mini1; min(Tcrt1(:,i))];
    maxi1 = [maxi1; max(Tcrt1(:,i))];
    avg_data1(i,1) = 0.5*(min(Tcrt1(:,i)) + max(Tcrt1(:,i)));
    std_PC1 = [std_PC1; std(Tcrt1(:,i))];
 
end

for i=1:size(Tcrt2,2)
    mini2 = [mini2; min(Tcrt2(:,i))];
    maxi2 = [maxi2; max(Tcrt2(:,i))];
    avg_data2(i,1) = 0.5*(min(Tcrt2(:,i)) + max(Tcrt2(:,i)));
    std_PC2 = [std_PC2; std(Tcrt2(:,i))];

end

for i=1:size(Tcrt3,2)
    mini3 = [mini3; min(Tcrt3(:,i))];
    maxi3 = [maxi3; max(Tcrt3(:,i))];
    avg_data3(i,1) = 0.5*(min(Tcrt3(:,i)) + max(Tcrt3(:,i)));
    std_PC3 = [std_PC3; std(Tcrt3(:,i))];

end

for i=1:size(Tcrt005,2)
    mini005 = [mini005; min(Tcrt005(:,i))];
    maxi005 = [maxi005; max(Tcrt005(:,i))];
    avg_data005(i,1) = 0.5*(min(Tcrt005(:,i)) + max(Tcrt005(:,i)));
    std_PC005 = [std_PC005; std(Tcrt005(:,i))];

end

% z_score = input('Specify the Z Score: ');
z_score = [1.96;2.58]; % 95,99 CI

for i = 1:length(z_score)
    CI_up(:,i) = avg_data + z_score(i)*(std_PC./sqrt(size(Tcrt,1)));
    CI_low(:,i) = avg_data - z_score(i)*(std_PC./sqrt(size(Tcrt,1)));
end


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
plot([1:size(Tcrt,2)], std_PC, 'bo','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',8); hold on
plot([1:size(Tcrt1,2)], std_PC1, 'rv','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',8); hold on
plot([1:size(Tcrt2,2)], std_PC2, 'gp','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',8); hold on
plot([1:size(Tcrt3,2)], std_PC3, 'md','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1],'MarkerSize',8); hold on
xlabel('Principal Components')
ylabel('Standard Deviation')
legend ('Normal Operating Data','Fault # 1 Operating Data','Fault # 2 Operating Data',...
    'Fault # 3 Operating Data');
grid on
subplot(2,1,2);
plot([1:size(Tcrt,2)], std_PC, 'bo','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 1],'MarkerSize',8); hold on
plot([1:size(Tcrt1,2)], std_PC1, 'rv','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 0],'MarkerSize',8); hold on
plot([1:size(Tcrt2,2)], std_PC2, 'gp','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[0 0 0],'MarkerSize',8); hold on
plot([1:size(Tcrt3,2)], std_PC3, 'md','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor',[1 0 1],'MarkerSize',8); hold on
grid on
axis([1 10 1 3])
xlabel('Principal Components')
ylabel('Standard Deviation')
legend ('Normal Operating Data','Fault # 1 Operating Data','Fault # 2 Operating Data',...
    'Fault # 3 Operating Data');

hold off;

%% Area/Centroid of Trapezoid
%Center of area or center of mass for a uniform lamina
centroid = []; centroid1 = [];centroid2 = [];centroid3 = [];
area = []; area1 = [];area2 = [];area3 = [];
for i=1:size(Tcrt,2)-1
    bb = max(Tcrt(:,i)) - min(Tcrt(:,i));
    aa = max(Tcrt(:,i + 1)) - min(Tcrt(:,i + 1));
    centroid(i) = (bb + 2*aa)/(bb + aa)/3;
    area(i) = (bb + aa)/2;
end
for i=1:size(Tcrt1,2)-1
    bb = max(Tcrt1(:,i)) - min(Tcrt1(:,i));
    aa = max(Tcrt1(:,i + 1)) - min(Tcrt1(:,i + 1));
    centroid1(i) = (bb + 2*aa)/(bb + aa)/3;
    area1(i) = (bb + aa)/2;
end
for i=1:size(Tcrt2,2)-1
    bb = max(Tcrt2(:,i)) - min(Tcrt2(:,i));
    aa = max(Tcrt2(:,i + 1)) - min(Tcrt2(:,i + 1));
    centroid2(i) = (bb + 2*aa)/(bb + aa)/3;
    area2(i) = (bb + aa)/2;
end
for i=1:size(Tcrt3,2)-1
    bb = max(Tcrt3(:,i)) - min(Tcrt3(:,i));
    aa = max(Tcrt3(:,i + 1)) - min(Tcrt3(:,i + 1));
    centroid3(i) = (bb + 2*aa)/(bb + aa)/3;
    area3(i) = (bb + aa)/2;
end
% Centroid plotting 
figure()
subplot(2,1,1);
aa = (1:size(Tcrt,2))-0.5;    
centroidaxis = aa(1,2:size(aa,2)); % used 2:size..to omit first x axis label
plot(centroidaxis,centroid,'b*');hold on
plot(centroidaxis,area,'r*');hold on

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

% Area plotting 

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
    
    figure
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
figure(2)
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

%% Plotting normal data 000 (1 std dev) or 005 (0.5 std dev) on normal data 00

figure()
subplot(2,1,1);
plot([1:size(Tcrt,2)], avg_data, 'bo-','linewidth',1.15); hold on
plot([1:size(Tcrt1,2)], avg_data1, 'r*-','linewidth',1.15);hold on
plot([1:size(Tcrt2,2)], avg_data2, 'g>-','linewidth',1.15);hold on
plot([1:size(Tcrt3,2)], avg_data3, 'md-','linewidth',1.15);hold on;
plot([1:size(Tcrt005,2)], avg_data005, 'bd-','linewidth',1.15);hold on;

plot([1:size(Tcrt,2)], CI_up(:,2),'k-.','linewidth',1);%99 CI
plot([1:size(Tcrt,2)], CI_low(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_up(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_low(:,2),'k-.','linewidth',1);%99 CI
xlabel('Principal Components')
ylabel('Mean and CI')
legend ('Normal Operating Data','Fault # 1','Fault # 2',...
    'Fault # 3','NO 0.5 SD','99% CI','95% CI');
grid on
%figure
subplot(2,1,2);
plot([1:size(Tcrt,2)], avg_data, 'bo-','linewidth',1.15); hold on
plot([1:size(Tcrt1,2)], avg_data1, 'r*-','linewidth',1.15);hold on
plot([1:size(Tcrt2,2)], avg_data2, 'g>-','linewidth',1.15);hold on
plot([1:size(Tcrt3,2)], avg_data3, 'md-','linewidth',1.15);hold on;
plot([1:size(Tcrt005,2)], avg_data005, 'bd-','linewidth',1.15);hold on;

plot([1:size(Tcrt,2)], CI_up(:,2),'k-.','linewidth',1);%99 CI
plot([1:size(Tcrt,2)], CI_low(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_up(:,1),'r-.','linewidth',1);%95 CI
plot([1:size(Tcrt,2)], CI_low(:,2),'k-.','linewidth',1);%99 CI
grid on
axis([1 10 -0.7 1.2])
xlabel('Principal Components')
ylabel('Mean and CI')
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

%%
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
%%