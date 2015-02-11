%% Plot for variables with most variation 

figure 
% for i=1:size(Xcrt,2)
% % subplot(4,1,1)
% 
% plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,i))+1.96*std(Xcrt(:,i)),1,size(Xcrt,1)),'k');hold on
% plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,i))-1.96*std(Xcrt(:,i)),1,size(Xcrt,1)),'k');hold on
% 
% plot(Xcrt(:,i),'bo'); hold on
% plot(Xcrt3(:,i),'ro'); hold on
% line([160 160], [-3 3],'Color',[1 0 1]); hold on
% axis([1 960 -3 3]);
% temp = strcat('variable #','-',int2str(i));
% xlabel(temp)
% 
% pause(5)
% hold off
% end

% plot most varying variables in normalized data  in fault # 3
subplot(5,2,1)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,7))+1.96*std(Xcrt(:,7)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,7))-1.96*std(Xcrt(:,7)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,7),'bo'); hold on
plot(Xcrt3(:,7),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('Reactor Pressure')

subplot(5,2,2)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,11))+1.96*std(Xcrt(:,11)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,11))-1.96*std(Xcrt(:,11)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,11),'bo'); hold on
plot(Xcrt3(:,11),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('Product Separator Temperature')

subplot(5,2,3)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,13))+1.96*std(Xcrt(:,13)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,13))-1.96*std(Xcrt(:,13)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,13),'bo'); hold on
plot(Xcrt3(:,13),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('Product Separator Pressure')

subplot(5,2,4)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,16))+1.96*std(Xcrt(:,16)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,16))-1.96*std(Xcrt(:,16)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,16),'bo'); hold on
plot(Xcrt3(:,16),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('Stripper Pressure')

subplot(5,2,5)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,18))+1.96*std(Xcrt(:,18)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,18))-1.96*std(Xcrt(:,18)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,18),'bo'); hold on
plot(Xcrt3(:,18),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('Stripper Temperature')

subplot(5,2,6)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,19))+1.96*std(Xcrt(:,19)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,19))-1.96*std(Xcrt(:,19)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,19),'bo'); hold on
plot(Xcrt3(:,19),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('Stripper steam flow')

subplot(5,2,7)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,20))+1.96*std(Xcrt(:,20)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,20))-1.96*std(Xcrt(:,20)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,20),'bo'); hold on
plot(Xcrt3(:,20),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('Compressor Work')


subplot(5,2,8)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,24))+1.96*std(Xcrt(:,24)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,24))-1.96*std(Xcrt(:,24)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,24),'bo'); hold on
plot(Xcrt3(:,24),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('E feed flow valve')


subplot(5,2,9)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,27))+1.96*std(Xcrt(:,27)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,27))-1.96*std(Xcrt(:,27)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,27),'bo'); hold on
plot(Xcrt3(:,27),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('Compressor recycle valve')


subplot(5,2,10)

plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,31))+1.96*std(Xcrt(:,31)),1,size(Xcrt,1)),'k');hold on
plot(1:size(Xcrt,1),repmat(mean(Xcrt(:,31))-1.96*std(Xcrt(:,31)),1,size(Xcrt,1)),'k');hold on
plot(Xcrt(:,31),'bo'); hold on
plot(Xcrt3(:,31),'ro'); hold on
line([160 160], [-3 3],'Color',[1 0 1]); hold on
axis([1 960 -3 3]);
xlabel('Stripper Steam valve')

%% CUMSUM
% cumsum(A,1) works along the rows of A and returns the cumulative sum of each column.
K = 2;  
c_plus = zeros(size(Tcrt,1),size(Tcrt,2));
c_minus = zeros(size(Tcrt,1),size(Tcrt,2));
avg_c_plus = zeros(size(Tcrt,1),size(Tcrt,2));
avg_c_minus = zeros(size(Tcrt,1),size(Tcrt,2));
for j=1:size(Tcrt,2) 
    for i=2:size(Tcrt,1) %skipping the first row coz of i-1 inital condition in next line
        c_plus(i,j)= max(0,(Tcrt(i,j)-K+c_plus(i-1,j)));
        c_minus(i,j)= max(0,-K-(Tcrt(i,j)+c_minus(i-1,j)));
        avg_c_plus(i,j)= (0.5 * (abs(c_plus(i-1,j))+abs(c_plus(i,j))));
        avg_c_minus(i,j)= (0.5 * (abs(c_minus(i-1,j))+abs(c_minus(i,j))));
    end
end

c_plus1 = zeros(size(Tcrt,1),size(Tcrt,2));
c_minus1 = zeros(size(Tcrt,1),size(Tcrt,2));
for j=1:size(Tcrt,2) 
    for i=2:size(Tcrt,1) %skipping the first row coz of i-1 inital condition in next line
        c_plus1(i,j)= max(0,(Tcrt1(i,j)-K+c_plus1(i-1,j)));
        c_minus1(i,j)= max(0,(-K-Tcrt1(i,j)+c_minus1(i-1,j)));
    end
end

c_plus3 = zeros(size(Tcrt,1),size(Tcrt,2));
c_minus3 = zeros(size(Tcrt,1),size(Tcrt,2));
for j=1:size(Tcrt,2) 
    for i=2:size(Tcrt,1) %skipping the first row coz of i-1 inital condition in next line
        c_plus3(i,j)= max(0,(Tcrt3(i,j)-K+c_plus3(i-1,j)));
        c_minus3(i,j)= max(0,-K-(Tcrt3(i,j)+c_minus3(i-1,j)));
    end
end
cusum = cumsum(Tcrt,1);
figure;plot(1:size(cusum,1),cusum(:,1),'go');hold on
plot(1:size(cusum,1),mean(Tcrt(:,1))+1.96*std(Tcrt(:,1)),'k-');
plot(1:size(cusum,1),mean(Tcrt(:,1))-1.96*std(Tcrt(:,1)),'k-');

cusum1 = cumsum(Tcrt,1);
figure;plot(cusum1(:,1));

cusum1 = cumsum(Xcrt,1);
figure;plot(cusum1(:,1));

cusum3 = cumsum(Tcrt3,1);
figure;plot(cusum3(:,1));

% envelope of NO process over the faulty data 
figure;plot(1:size(c_plus,1),c_plus(:,1),'k-');hold on
plot(avg_c_plus(:,1),'b'); hold on
plot(1:size(c_plus,1),c_minus(:,1),'r-'); hold on
plot(avg_c_minus(:,1),'b');

figure;plot(1:size(c_plus,1),c_plus1(:,1),'k-');hold on
plot(avg_c_plus(:,1),'b'); hold on
plot(1:size(c_plus,1),c_minus1(:,1),'r-'); hold on
plot(avg_c_minus(:,1),'b');

figure;plot(1:size(c_plus,1),c_plus3(:,1),'k-');hold on
plot(avg_c_plus(:,1),'b'); hold on
plot(1:size(c_plus,1),c_minus3(:,1),'r-'); hold on
plot(avg_c_minus(:,1),'b');

% how cumsum is different from actual scores
figure
plot(Tcrt(:,1),'g'); hold on 
plot(cumsum(Tcrt(:,1))*0.1,'r')

%Now since the CUSUM detection is delayed as seen in plots above, I want to
%see how the plots vary with absolute CUSUM values 
figure
for i=1:4
    
plot(cumsum(abs(Tcrt(:,i)))*0.1,'b');hold on
plot(cumsum(abs(Tcrt3(:,i)))*0.1,'g');hold on
plot(cumsum(abs(Tcrt1(:,i)))*0.1,'r');hold off
title(i);
pause(20)
end

% for normalized data 
figure
for i=1:4
    
plot(cumsum((Xcrt(1:200,i)))*0.1,'b');hold on
plot(cumsum((Xcrt3(1:200,i)))*0.1,'g');hold on
plot(cumsum((Xcrt1(1:200,i)))*0.1,'r');hold off
title(i);
pause(30)
end

%% EWMA

lamda = 0.2;  L = 2.962; 
ewma = zeros(size(Tcrt,1),size(Tcrt,2));
ewma_up = zeros(size(Tcrt,1),size(Tcrt,2));
ewma_low = zeros(size(Tcrt,1),size(Tcrt,2));
for j=1:size(Tcrt,2) 
    for i=2:size(Tcrt,1) %skipping the first row coz of i-1 inital condition in next line
        ewma(i,j) = lamda*Tcrt(i,j)+(1-lamda)*ewma(i-1,j); 
        ewma_up(i,j) = L*2.42*((lamda/(2-lamda))*sqrt(1-(1-lamda))); 
        ewma_low(i,j) = -L*2.42*((lamda/(2-lamda))*sqrt(1-(1-lamda))); 
    end
end

plot(1:size(ewma,1),ewma(:,1));hold on
plot(1:size(ewma,1),ewma_up(:,1));hold on
plot(1:size(ewma,1),ewma_low(:,1));

%% Scree plots 

k=1; % fault #
for i=1:1 % this is the PC#
figure()

temp = []; temp1 = [];temp2 = [];temp3 = [];

eval(['temp','=','Tcrt',num2str(k),'(:,','i)',';']);
eval(['temp1','=','Tcrt',num2str(k),'(:,','i+1)',';']);
eval(['temp2','=','Pref','(:,','i)',';']);
eval(['temp3','=','Pref','(:,','i+1)',';']);

% plot(temp,'ro','MarkerSize',2); hold on 
%if you dont want sample # on plot uncomment the code above

% code below for adding sample # on plot
j=(1:960)'; % this is the sample size used for plotting 
% strValues1 = strtrim(int2str(j));
% text(1:size(Tcrt,1),temp,strValues1,'VerticalAlignment','bottom');hold on
% text(1:size(Tcrt,1),temp1,strValues1,'VerticalAlignment','bottom');hold on
plot3(j,temp,temp1,'ro','MarkerFaceColor','r','Markersize',2);hold on
xlabel('sample #');
ylabel('Scores PC 1');
zlabel('Scores PC 2');

figure()
j=(1:33)';
plot3(j,temp2,temp3,'ro','MarkerFaceColor','r','Markersize',2);hold on
xlabel('sample #');
ylabel('Loading on PC 1');
zlabel('Loading on PC 2');
% temp = strcat('PC #','-',int2str(i));
% title(temp);
% axis([1 10 -6 6]);
% set(gca,'xgrid','on')
% pause(1)
hold off
end

%% Row by row std dev of PCA scores

std_norm = [];
for i = 1:size(Tcrt,1)
    std_norm = [std_norm; std(Tcrt(i,:))];
end

std_norm1 = [];
for i = 1:size(Tcrt1,1)
    std_norm1 = [std_norm1; std(Tcrt1(i,:))];
end
std_norm3 = [];
for i = 1:size(Tcrt3,1)
    std_norm3 = [std_norm3; std(Tcrt3(i,:))];
end

figure()
plot(std_norm,'bo'); hold on
plot(std_norm3,'ro'); hold on
plot(std_norm1,'go'); hold on
line([160 160], [0 5],'Color',[1 0 1]); hold on
plot(1:size(Tcrt,1),mean(std_norm)+1.96*std(std_norm),'k');hold on
plot(1:size(Tcrt,1),mean(std_norm)-1.96*std(std_norm),'k');hold on


figure()
plot(std_norm,'bo'); hold on
plot(std_norm3,'ro'); hold on
line([160 160], [0 4],'Color',[1 0 1]); hold on
plot(1:size(Tcrt,1),mean(std_norm)+1.96*std(std_norm),'k');hold on
plot(1:size(Tcrt,1),mean(std_norm)-1.96*std(std_norm),'k');hold on
hold off

%% Row by Row for all faults with pause  
%%find std_norm for all faults 
    for i=1:21 
        eval(['std_norm',num2str(i),'=','[]'';']);
        for j=1:size(Tcrt,1)
        eval(['std_norm',num2str(i),'=','[','std_norm',num2str(i),';',...
            'std(Tcrt',num2str(i),'(j,:))]',';']);
        end
    end
  
    
for i=3:3 % this is the fault#
figure()
plot(std_norm,'bo','MarkerFaceColor','b','Markersize',2); hold on

temp = [];
eval(['temp','=','std_norm',num2str(i),';']);
% plot(temp,'ro','MarkerSize',2); hold on 
%if you dont want sample # on plot uncomment the code above

% code below for adding sample # on plot
j=(1:960)'; % this is the sample size used for plotting 
strValues1 = strtrim(int2str(j));
text(1:size(Tcrt,1),temp,strValues1,'VerticalAlignment','bottom');
hold on;
plot(temp,'ro','MarkerFaceColor','r','Markersize',2);
hold on;

plot(1:size(Tcrt,1),mean(std_norm)+1.96*std(std_norm),'g');hold on
plot(1:size(Tcrt,1),mean(std_norm)-1.96*std(std_norm),'g');hold on
line([160 160], [0 2.5],'Color',[1 0 1]); hold on
temp = strcat('Fault #','-',int2str(i));
title(temp);
% axis([1 10 -6 6]);
% set(gca,'xgrid','on')
% pause(1)
hold off
end



%% std dev of each PC, this is for just a particular fault

for i=1:10 % This is the # of PC you want to check for 
plot(1:size(Tcrt,1),Tcrt3(:,i),'bo','Markersize',2);hold on    %Change the Tcrt only here for fault
plot(1:size(Tcrt,1),1.96*std(Tcrt(:,i)),'r-');hold on
plot(1:size(Tcrt,1),-1.96*std(Tcrt(:,i)),'r-');hold on
plot([160 160],[1.96*std(Tcrt(:,i)) -1.96*std(Tcrt(:,i))],'m'); hold on % how to add vertical line? 
temp = strcat('Principal Component','-',int2str(i));
xlabel(temp);
pause(1)
hold off
end

%% Parallel Coordinates - std dev of each PC,this is for just a particular fault
limit_up= std(Tcrt(:,:));
for i=220:220 % This is the # of time samples you want to check for 

plot(1:size(Tcrt,2),Tcrt3(i,:),'b*','MarkerSize',10);hold on    %Change the Tcrt only here for fault
plot(1:size(Tcrt,2),1.96*limit_up,'r');hold on
plot(1:size(Tcrt,2),-1.96*limit_up,'r');hold on
% plot([160 160],[1.96*std(Tcrt(:,i)) -1.96*std(Tcrt(:,i))],'m'); hold on % how to add vertical line? 

temp = strcat('sample #','-',int2str(i));
title(temp);
axis([1 10 -6 6]);
set(gca,'xgrid','on')
pause(4)
hold off
end

%%
figure()
subplot(4,1,1);
plot(1:size(Tcrt,2),Tcrt(:,:),'b*','MarkerSize',3);hold on    %Change the Tcrt only here for fault
plot(1:size(Tcrt,2),1.96*limit_up,'r');hold on
plot(1:size(Tcrt,2),-1.96*limit_up,'r');hold on
title('Normal');
subplot(4,1,2);
plot(1:size(Tcrt,2),Tcrt3(:,:),'b*','MarkerSize',3);hold on    %Change the Tcrt only here for fault
plot(1:size(Tcrt,2),1.96*limit_up,'r');hold on
plot(1:size(Tcrt,2),-1.96*limit_up,'r');hold on
title('Fault 3');
subplot(4,1,3);
plot(1:size(Tcrt,2),Tcrt9(:,:),'b*','MarkerSize',3);hold on    %Change the Tcrt only here for fault
plot(1:size(Tcrt,2),1.96*limit_up,'r');hold on
plot(1:size(Tcrt,2),-1.96*limit_up,'r');hold on
title('Fault 9');
subplot(4,1,4);
plot(1:size(Tcrt,2),Tcrt1(:,:),'b*','MarkerSize',3);hold on    %Change the Tcrt only here for fault
plot(1:size(Tcrt,2),1.96*limit_up,'r');hold on
plot(1:size(Tcrt,2),-1.96*limit_up,'r');hold on
title('Fault 1');

%% Fault detection 
            
cutoff_up= (1.96*limit_up); % store upper 99% values in cutoff_up
cutoff_low= (-1.96*limit_up); % store lower 99% values in cutoff_low
sum_norm=[]; sum_05sd=[];sum_1sd=[];sum_2sd=[];

faultsmatrix = zeros(size(Tcrt,1),size(Tcrt,2));   
faultsmatrix_norm = zeros(size(Tcrt,1),size(Tcrt,2));
faultsmatrix_05sd = zeros(size(Tcrt,1),size(Tcrt,2));
faultsmatrix_1sd = zeros(size(Tcrt,1),size(Tcrt,2));
faultsmatrix_2sd = zeros(size(Tcrt,1),size(Tcrt,2));

    for i=1:21                          % # of fault data sets
        for j=1:size(Tcrt,2)            % # of columns in Tcrt
            eval(['temp','=','Tcrt',num2str(i),';']);
            for k=1:size(Tcrt,1)
                if (temp(k,j)> cutoff_up(1,j)) || (temp(k,j)< cutoff_low(1,j))
                    faultsmatrix(k,j)=1; 
                end    
            end
        end
        eval(['faultsmatrix',num2str(i),'=','faultsmatrix',';']);
        eval(['sum_faults','(i,:)','=','sum(faultsmatrix)',';']); %column # is the fault
        faultsmatrix = zeros(size(Tcrt,1),size(Tcrt,2));
        temp = zeros(size(Tcrt,1),size(Tcrt,2));
    end
% row of sum_faults is the sum of # of faults detected 
% 

for i=1:21   
    eval(['temp','=','faultsmatrix',num2str(i),';']);    
    for j=1:size(faultsmatrix,2) %# of columns 
        for k=2:size(Tcrt,1) % # of rows, start from 2 because j-1 in the next step
            if (temp(k,j)> 0) && (temp(k,j)== temp(k-1,j))
                sample(k,j) = k; %sample(k-1,j) = k-1, k-1 is the first fault
            end
        end
    end
    eval(['sample',num2str(i),'=','sample',';'])
    sample = zeros(size(Tcrt,1),size(Tcrt,2));    
end   
    

% %   for i=1:size(cutoff_up,1) % all rows
% %             if (avg_data(i,1)> cutoff_up(i,1)) || (avg_data(i,1)< cutoff_low(i,1))
% %             faultsmatrix_norm(i,1)=1;
% %             end       
% %             sum_norm = sum(faultsmatrix_norm);
% %   end
% %   for i=1:size(cutoff_up,1) % all rows
% %             if (avg_data05sd(i,1)> cutoff_up(i,1)) || (avg_data05sd(i,1)< cutoff_low(i,1))
% %             faultsmatrix_05sd(i,1)=1;
% %             end
% %             sum_05sd = sum(faultsmatrix_05sd);    
% %   end
% %   for i=1:size(cutoff_up,1) % all rows
% %             if (avg_data1sd(i,1)> cutoff_up(i,1)) || (avg_data1sd(i,1)< cutoff_low(i,1))
% %             faultsmatrix_1sd(i,1)=1;
% %             end
% %             sum_1sd = sum(faultsmatrix_1sd);                 
% %   end
% %   for i=1:size(cutoff_up,1) % all rows
% %             if (avg_data2sd(i,1)> cutoff_up(i,1)) || (avg_data2sd(i,1)< cutoff_low(i,1))
% %             faultsmatrix_2sd(i,1)=1;
% %             end
% %             sum_2sd = sum(faultsmatrix_2sd);                   
% %   end

%% Contribution Plots

for i=168:180 % This is the # of time samples you want to check for 

temp = strcat('sample #',int2str(i));

% contribution plots for Hotelling's T2
ConT1 = (Xcrt3(i, :)*Pref(:, 1:a)*sqrt(inv(diag(Eref(1:a))))*Pref(:, 1:a)').^2;
subplot(2,2,1);
plot(1:size(ConT1,2),T2c,'k--');
bar(1:size(ConT1,2),ConT1,'r');

xlabel('Process Variables'); ylabel('Contribution to T2');
legend(temp);
axis([0 34 -1 5])
grid minor
grid on

% contribution plots for Q
ConQ1 = (Xcrt3(i,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')).^2;
subplot(2,2,2);
plot(1:size(ConT1,2),Qc,'k--');
bar(ConQ1,'r');xlabel('Process Variables'); ylabel('Contribution to Q');

legend(temp);
axis([0 34 -1 8])
grid minor
grid on

subplot(2,2,[3,4])
plot(1:size(Tcrt,2),Tcrt3(i,:),'b*','MarkerSize',10);hold on    %Change the Tcrt only here for fault
plot(1:size(Tcrt,2),1.96*limit_up,'r');hold on
plot(1:size(Tcrt,2),-1.96*limit_up,'r');hold on
% plot([160 160],[1.96*std(Tcrt(:,i)) -1.96*std(Tcrt(:,i))],'m'); hold on % how to add vertical line? 

temp = strcat('sample #','-',int2str(i));
title(temp);
axis([1 10 -6 6]);
set(gca,'xgrid','on')
pause(1)
hold off
end

ConT1 = (Xcrt3(:, :)*Pref(:, 1:a)*sqrt(inv(diag(Eref(1:a))))*Pref(:, 1:a)').^2;
ConQ1 = (Xcrt3(:,:)*(eye(p)-Pref(:,1:a)*Pref(:,1:a)')).^2;
temp=sum(ConT1,2);
temp1=sum(ConQ1,2);
figure
plot(temp)
figure
plot(temp1)

figure
plot(1:size(ConT1,2),ConT1(100:180,:),'g');hold on 
plot(1:size(ConT1,2),ConT1(170:178,:),'r');hold on 
figure
plot(1:size(ConT1,2),ConQ1(100:180,:),'g');hold on 
plot(1:size(ConT1,2),ConQ1(170:178,:),'r');hold on 

%%
figure
subplot(4,2,1)
bar(ConT1(160,:))
title('T2');
xlabel('sample#160');
axis([1 33 0 2])
subplot(4,2,3)
bar(ConT1(170,:))
xlabel('sample#170');
axis([1 33 0 2])
subplot(4,2,5)
bar(ConT1(178,:))
xlabel('sample#178');
axis([1 33 0 2])
subplot(4,2,7)
bar(ConT1(171,:))
xlabel('sample#171');
axis([1 33 0 2])

subplot(4,2,2)
bar(ConQ1(160,:))
xlabel('sample#160');
ylabel('Q');
axis([1 33 0 2])
subplot(4,2,4)
bar(ConQ1(170,:))
xlabel('sample#170');
axis([1 33 0 2])
subplot(4,2,6)
bar(ConQ1(178,:))
xlabel('sample#178');
axis([1 33 0 2])
subplot(4,2,8)
bar(ConQ1(171,:))
xlabel('sample#171');
axis([1 33 0 2])
%% plot T2 and Q for the faulty samples 

nzero = sort((nonzeros(sample3(:,1))),'ascend');

figure 
for i=1:size(nzero,1)%160:180 %1:size(nzero,1)
subplot(2,1,1)
j=nzero(i,1);

plot(1:size(ConT1,2),ConT1(j,:),'r'); hold on
plot(1:size(ConT1,2),ConT1(100:160,:),'g'); hold on
% plot(1:size(ConT1,2),ConT1(i,:),'b');
set(gca,'XTick',[1:33]);
grid minor 
hold off 

subplot(2,1,2)

plot(1:size(ConT1,2),ConQ1(j,:),'r'); hold on
plot(1:size(ConT1,2),ConQ1(100:160,:),'g'); hold on
% plot(1:size(ConT1,2),ConQ1(i,:),'b');
set(gca,'XTick',[1:33]);
grid minor 

temp = strcat('sample #','-',int2str(j));
title(temp);
pause(1)
hold off
end

%%
%for faulty points
figure
for i=1:size(nzero,1)
j=nzero(i,1);
bar(ConT1(j,:),'r')
temp = strcat('sample #','-',int2str(j));
axis([1 33 0 2])
title(temp);
pause(1)
end

% for normal samples
figure
for i=65:100
bar(ConQ1(i,:),'r')
temp = strcat('sample #','-',int2str(i));
axis([1 33 0 2])
title(temp);
pause(1)
end

% subplot(2,1,2)
% plot(1:size(ConT1,2),ConQ1(:,:))
%%
figure
subplot(2,2,1)
plot(d00_te(:,11),'bo'); hold on
% subplot(2,2,2)
plot(d03_te(:,11),'ro')

%% Old Fault detection step method 
ss=[];kk = 1;yy= 1; %kk was 30, yy is the sample size
max_sample = 160; sss=max_sample/kk; %max_sample was 960
fault_det = zeros(sss,size(Tcrt,2));
fault_det05sd = zeros(sss,size(Tcrt,2));  
fault_det1sd = zeros(sss,size(Tcrt,2)); 
fault_det2sd = zeros(sss,size(Tcrt,2)); 
fault_det3 = zeros(sss,size(Tcrt,2));  
fault_det9 = zeros(sss,size(Tcrt,2)); 
cutoff_up= CI_up(:,2)'; % store upper 99% values in cutoff_up
cutoff_low=CI_low(:,2)'; % store lower 99% values in cutoff_low

%to get sample means 
    for j=1:sss
        fault_det(j,:) = mean(Tcrt(kk:kk+yy,:)); %each row is mean for given step
        fault_det3(j,:) = mean(Tcrt3(kk:kk+yy,:));
        fault_det05sd(j,:) = mean(Tcrt05sd(kk:kk+yy,:));
        fault_det1sd(j,:) = mean(Tcrt1sd(kk:kk+yy,:));
        fault_det2sd(j,:) = mean(Tcrt2sd(kk:kk+yy,:));
        fault_det9(j,:) = mean(Tcrt9(kk:kk+yy,:));
        ss(j,1)= kk;
       

        if kk < 160 %max_sample-kk should be 930
        kk= kk + yy ; 
        end
    end
    
sum_det= [];
fault_det_matrix = zeros(size(fault_det,1),size(fault_det,2));
for j=1:size(fault_det,1) % all rows
     for i=1:size(cutoff_up,2) % all columns
            if (fault_det(j,i)> cutoff_up(1,i)) || (fault_det(j,i)< cutoff_low(1,i))
            fault_det_matrix(j,i)=1;
            end       
            sum_det = sum(fault_det_matrix,2);
     end 
end

fault_det_matrix = zeros(size(fault_det,1),size(fault_det,2));
for j=1:size(fault_det,1) % all rows
     for i=1:size(cutoff_up,2) % all columns
            if (fault_det3(j,i)> cutoff_up(1,i)) || (fault_det3(j,i)< cutoff_low(1,i))
            fault_det_matrix(j,i)=1;
            end       
            sum_det3 = sum(fault_det_matrix,2);
     end 
end
fault_det_matrix = zeros(size(fault_det,1),size(fault_det,2));
for j=1:size(fault_det,1) % all rows
     for i=1:size(cutoff_up,2) % all columns
            if (fault_det05sd(j,i)> cutoff_up(1,i)) || (fault_det05sd(j,i)< cutoff_low(1,i))
            fault_det_matrix(j,i)=1;
            end       
            sum_det05sd = sum(fault_det_matrix,2);
     end 
end
fault_det_matrix = zeros(size(fault_det,1),size(fault_det,2));
for j=1:size(fault_det,1) % all rows
     for i=1:size(cutoff_up,2) % all columns
            if (fault_det1sd(j,i)> cutoff_up(1,i)) || (fault_det1sd(j,i)< cutoff_low(1,i))
            fault_det_matrix(j,i)=1;
            end       
            sum_det1sd = sum(fault_det_matrix,2);
     end 
end
fault_det_matrix = zeros(size(fault_det,1),size(fault_det,2));
for j=1:size(fault_det,1) % all rows
     for i=1:size(cutoff_up,2) % all columns
            if (fault_det2sd(j,i)> cutoff_up(1,i)) || (fault_det2sd(j,i)< cutoff_low(1,i))
            fault_det_matrix(j,i)=1;
            end       
            sum_det2sd = sum(fault_det_matrix,2);
     end 
end
fault_det_matrix = zeros(size(fault_det,1),size(fault_det,2));

for j=1:size(fault_det,1) % all rows
     for i=1:size(cutoff_up,2) % all columns
            if (fault_det9(j,i)> cutoff_up(1,i)) || (fault_det9(j,i)< cutoff_low(1,i))
            fault_det_matrix(j,i)=1;
            end       
            sum_det9 = sum(fault_det_matrix,2);
     end 
end
