%% Fault detection entire data set
            
cutoff_up= CI_up(:,2); % store upper 99% values in cutoff_up
cutoff_low=CI_low(:,2); % store lower 99% values in cutoff_low
sum_norm=[]; sum_05sd=[];sum_1sd=[];sum_2sd=[];
faultsmatrix = zeros(size(cutoff_up,1),size(cutoff_up,2));   
faultsmatrix_norm = zeros(size(cutoff_up,1),size(cutoff_up,2));
faultsmatrix_05sd = zeros(size(cutoff_up,1),size(cutoff_up,2));
faultsmatrix_1sd = zeros(size(cutoff_up,1),size(cutoff_up,2));
faultsmatrix_2sd = zeros(size(cutoff_up,1),size(cutoff_up,2));

    for i=1:21                          % # of fault data sets
        for j=1:size(cutoff_up,1)            % # of columns in Tcrt
        eval(['temp_avg','=','avg_data',num2str(i),';']);
            if (temp_avg(j,1)> cutoff_up(j,1)) || (temp_avg(j,1)< cutoff_low(j,1))
            faultsmatrix(j,1)=1;
            sum_faults(:,i) = sum(faultsmatrix); %column # is the fault
            end    
        end
        faultsmatrix = zeros(size(cutoff_up,1),size(cutoff_up,2));  
    end
    
  for i=1:size(cutoff_up,1) % all rows
            if (avg_data(i,1)> cutoff_up(i,1)) || (avg_data(i,1)< cutoff_low(i,1))
            faultsmatrix_norm(i,1)=1;
            end       
            sum_norm = sum(faultsmatrix_norm);
  end
  for i=1:size(cutoff_up,1) % all rows
            if (avg_data05sd(i,1)> cutoff_up(i,1)) || (avg_data05sd(i,1)< cutoff_low(i,1))
            faultsmatrix_05sd(i,1)=1;
            end
            sum_05sd = sum(faultsmatrix_05sd);    
  end
  for i=1:size(cutoff_up,1) % all rows
            if (avg_data1sd(i,1)> cutoff_up(i,1)) || (avg_data1sd(i,1)< cutoff_low(i,1))
            faultsmatrix_1sd(i,1)=1;
            end
            sum_1sd = sum(faultsmatrix_1sd);                 
  end
  for i=1:size(cutoff_up,1) % all rows
            if (avg_data2sd(i,1)> cutoff_up(i,1)) || (avg_data2sd(i,1)< cutoff_low(i,1))
            faultsmatrix_2sd(i,1)=1;
            end
            sum_2sd = sum(faultsmatrix_2sd);                   
  end
  
%% Fault detection step method 
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

%% idea
TT = [ss,sum_det,sum_det05sd,sum_det1sd,sum_det2sd,sum_det3,sum_det9];
TT
size(TT)
fault_det_matrix = zeros(size(TT,1),size(TT,2));
sum_det_test= [];
for j=3:size(TT,2) % all columns
     for i=1:size(TT,1) % all rows
            if (TT(i,j)> TT(i,2)) 
            fault_det_matrix(i,j)=1;
            end       
            sum_det_test = sum(fault_det_matrix,1);
     end 
end
fault_det_matrix
LL = [ss,fault_det_matrix(:,3:size(TT,2))]
sum_det_test

