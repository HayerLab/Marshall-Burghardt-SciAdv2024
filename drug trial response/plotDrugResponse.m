% Script that plots FRET activity for the different drug conditions 
% Seph Marshall 200228

%% Control section 
clc; clear; 
control_root = 'F:example dataset';

datadir =('F:where to store data path');  

%first number is num sites per condition, 2nd is number of frames 
control_arr=zeros(20,60); 
 
k = 0; 
%specify the sites as determined by ND2 aquisition

for row = 1:4
    for col =1 
        for site = 1:5          
            k=k+1; 
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end 
    end 
end 


for i = 1:size(position,2)    
   imRatio_raw={}; 
    load([control_root,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
    
    for j = 1:60 %here j goes up to number of frames 
        
        temp_arr=imRatio_raw{1,j};
        
        control_arr(i,j)= nanmean(temp_arr,'all'); 
         
    end 
    
end




%% drug condition 1 section /
condition1_root = control_root; % change if different 

position = []; 
k = 0; 

condition1_arr=zeros(20, 60); 

for row = 1:4
    
    for col =1 
        for site =6:10
            k=k+1; 
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end 
    end 
end 




for i = 1:(size(position,2))
   
    load([condition1_root,filesep, 'data', filesep, position{i}, '_RatioData_raw.mat'],'imRatio_raw'); 
    
    
    for j = 1:60 %j goes to total number of frames
        
        temp_arr2=imRatio_raw{1,j};
        
        condition1_arr(i,j)= nanmean(temp_arr2,'all'); 
         
    end 
    
end 




%% normalize, get averages
% 

frames_pre_addition = 5; %number of frames pre drug addition, change as needed 
 
for i=1:size(control_arr,1)
     meanval = mean(control_arr(i,1:frames_pre_addition));
     control_arr(i,:)=control_arr(i,:)/meanval; 
 end 
 
 mean_control = mean(control_arr,1); 
 
 for i=1:size(control_arr,1)
    control_arr(i,:) = control_arr(i,:)./mean_control; 
 end 
 
stats_arrCNTRL = zeros(3,60); 

 for k = 1:60
     
     % fitting normal dist to each lag 
    pd = fitdist(control_arr(:,k),'Normal'); 
    stats_arrCNTRL(1,k) = pd.mu; 
    ci_CNTRL = paramci(pd); 
    
    %95 % confidence interval, upper and lower 
  stats_arrCNTRL(2,k) = ci_CNTRL(1,1); 
    stats_arrCNTRL(3,k)=ci_CNTRL(2,1); 
   
 end 



for i=1:size(condition1_arr,1)
    meanval = mean(condition1_arr(i,1:frames_pre_addition));
    condition1_arr(i,:)=condition1_arr(i,:)/meanval; 
end 


 
 for i=1:size(condition1_arr,1)
    condition1_arr(i,:) = condition1_arr(i,:)./mean_control; 
 end 
 
 
stats_arrTREAT = zeros(3,60); 

 for k = 1:60
     
     % fitting normal dist to each lag 
    pd = fitdist(condition1_arr(:,k),'Normal'); 
    stats_arrTREAT(1,k) = pd.mu; 
    ci_TREAT = paramci(pd); 
    
    %95 % confidence interval, upper and lower 
   
    stats_arrTREAT(2,k) = ci_TREAT(1,1); 
    stats_arrTREAT(3,k) = ci_TREAT(2,1); 
 end 

%% plot the data

f1=figure; 

xlim([0 60]); 


 xline(10,'--'); 

 axis square; 
 hold on; 
 
  title('POI, control vs condition 1');
 xlabel('TimePoint (1 min)'); 
 ylabel('Norm. POI'); 

for a=1:size(control_arr,1)
    plot([1:60],control_arr(a,:), 'Color','k','DisplayName','CNTRL'); 
 
end 


for c=1:size(condition1_arr,1)
    plot([1:60],condition1_arr(c,:), 'Color','r','DisplayName','+drug @5'); 
   
end 
 
hold off; 

f2= figure; 
 title('mean +/- CI'); 
 xlabel('TimePoint (1 min)'); 
 ylabel('Norm. POI'); 
 hold on; 
 ylim([0.9 1.25]);
 xlim([0 60]);
xline(5,'--'); 

axis square; 

 

plot(1:60,(stats_arrCNTRL(1,:)),'Color',[0,0,0], 'LineWidth', 3 ); 

plot(1:60,stats_arrCNTRL(2,:),'Color',[0,0,0] );
plot(1:60,stats_arrCNTRL(3,:),'Color',[0,0,0] );

plot(1:60,(stats_arrTREAT(1,:)),'Color',[1,0,0], 'LineWidth', 3 ); 
plot(1:60,(stats_arrTREAT(2,:)),'Color',[1,0,0] );
plot(1:60,(stats_arrTREAT(3,:)),'Color',[1,0,0] );


    save([datadir,filesep,'drug trial data.mat'],'control_arr', 'condition1_arr', 'stats_arrCNTRL', 'stats_arrTREAT'); 
    saveas(f1, ([datadir, filesep, 'graph.fig'])); 
     saveas(f1, ([datadir, filesep, 'graph.svg'])); 
    saveas(f2, ([datadir, filesep, 'graph_mean_CI.fig'])); 
    saveas(f2, ([datadir, filesep, 'graph_mean_CI.svg'])); 