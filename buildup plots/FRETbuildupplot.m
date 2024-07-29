% script that compiles multiple protrusion/retraction events inot one graph
% with mean and confidence intervals 

clear; 
clc; 

root = 'F:\example set 1\FRET buildup - retractions'; 
root2= 'F:\example set 2\FRET buildup - retractions'; 
root3 = 'F:\example set 3\FRET buildup - retractions'; 

datadir ='C:\save where you desire'; 

%first number is number of total events, second number is your frame range 
% ie -15:1:30 is 46 frames
velocities = zeros(30,46); 
FRET = zeros(30,46); 
myosin = zeros(30,46); 

%number of events from dataset 1
for i=1:10
    
   
    load([root, filesep, num2str(i), filesep, 'retraction_statistics.mat']); 
    
    
    velocities(i,:) = vel_arr; 
    
    % here specify edge depth you want to plot
    %FRET_temp has 6 depths stored, and you choose the one you desire 
    for w=1:size(FRET_temp{2,2},2)
        FRET_avg(1,w) =mean(FRET_temp{2,2}(:,w));     
    end 

    FRET(i,:) = FRET_avg; 
    
    %here specify which edge depth you wish to plot 
    for w=1:size(myosin_temp{2,5},2)
        myosin_avg(1,w) =mean(myosin_temp{2,5}(:,w));     
    end 

    myosin(i,:) = myosin_avg; 
    
    
end

%number of events from dataset 2 
for i=11:20
    
    load([root2, filesep, num2str(i-10), filesep, 'retraction_statistics.mat']); 
   
    velocities(i,:) = vel_arr; 
  
    % here specify edge depth you want to plot 
    for w=1:size(FRET_temp{2,2},2)
        FRET_avg(1,w) =mean(FRET_temp{2,2}(:,w));     
    end 

    FRET(i,:) = FRET_avg; 
 
    for w=1:size(myosin_temp{2,5},2)
        myosin_avg(1,w) =mean(myosin_temp{2,5}(:,w));     
    end 

    myosin(i,:) = myosin_avg; 
    
  end 
% %     
 
% number of events in dataset 3 
for i=21:30
    
    load([root3, filesep, num2str(i-20), filesep, 'retraction_statistics.mat']); 
    
    velocities(i,:) = vel_arr; 
    
    
    % here specify edge depth you want to plot 
    for w=1:size(FRET_temp{2,2},2)
        FRET_avg(1,w) =mean(FRET_temp{2,2}(:,w));     
    end 

    FRET(i,:) = FRET_avg; 
    
    for w=1:size(myosin_temp{2,5},2)
        myosin_avg(1,w) =mean(myosin_temp{2,5}(:,w));     
    end 

    myosin(i,:) = myosin_avg; 
    
    
end 



 for k = 1:46 % -15:30 lags, respectively 
     
     % fitting normal dist to each timepoint 
    pd = fitdist(velocities(:,k),'Normal'); 
    velocity_arr(1,k) = pd.mu; 
    ci = paramci(pd); 
    
    velocity_arr(2,k) = ci(1,1); 
    velocity_arr(3,k)=ci(2,1); 
    
     pd2 = fitdist(FRET(:,k),'Normal'); 
    FRET_arr(1,k) = pd2.mu; 
    ci2 = paramci(pd2); 
    
    %95 % confidence interval, upper and lower 
    FRET_arr(2,k) = ci2(1,1); 
    FRET_arr(3,k)=ci2(2,1); 
    
       pd3 = fitdist(myosin(:,k),'Normal'); 
     myosin_arr(1,k) = pd3.mu; 
     ci3 = paramci(pd3); 
%     
%     %95 % confidence interval, upper and lower 
    myosin_arr(2,k) = ci3(1,1); 
    myosin_arr(3,k)=ci3(2,1); 
 end 
% 
%% plotting figures 

x = -15:1:30; 

f=figure; 

hold on; 
yyaxis left; 
% commented out option to plot edge velocity  

% ylabel('Edge Velocity (um/min)','Color',[0,0,1]);
% ylim([-17 17]); 
% yticks([ -15.38 -12.82 -10.256 -7.692 -5.128 -2.564 0 2.564 5.128 7.692 10.256 12.82 15.38 ]); 
% yticklabels({ '-12' '-10', '-8', '-6', '-4', '-2', '0', '2', '4', '6', '8', '10' '12'}); 

% plot(x,velocity_arr(1,:),'-','Color',[0,0,1], 'LineWidth', 3 );
% plot(x,velocity_arr(2,:),'-','Color',[0,0,1], 'LineWidth', 1 );
% plot(x,velocity_arr(3,:),'-','Color',[0,0,1], 'LineWidth', 1 );

% xticks([-19.2 -16.8 -14.4 -12.0 -9.6 -7.2 -4.8 -2.4 0 2.4 4.8 7.2 9.6 12.0 14.4 16.8 19.2])
% xticklabels({'-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8'}); 
% 
xlabel('Time (min)'); 
 xline(0,'--'); 
 yline(0,'--k'); 

 xticks([ -14.4 -12.0 -9.6 -7.2 -4.8 -2.4 0 2.4 4.8 7.2 9.6 12.0 14.4 16.8 19.2 21.6 24.0 26.4 28.8 ])
 xticklabels({ '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8' '' '10' '' '12'}); 

ylabel(' Norm. MLC-mRuby3'); 
% ylim([0 0.8]);
 

plot(x,myosin_arr(1,:),'-','Color',[0,0,1], 'LineWidth', 3 );
plot(x,myosin_arr(2,:),'-','Color',[0,0,1], 'LineWidth', 1 );
plot(x,myosin_arr(3,:),'-','Color',[0,0,1], 'LineWidth', 1 );




yyaxis right; 
ylabel('DORA-RhoB'); 
ylim([0.6 1.4]);
title('retraction progression (RhoB 1.95 um, MLC 6.5 um');

plot(x,FRET_arr(1,:),'-','Color',[1,0,0], 'LineWidth', 3 );
plot(x,FRET_arr(2,:),'-','Color',[1,0,0], 'LineWidth', 1 );
plot(x,FRET_arr(3,:),'-','Color',[1,0,0], 'LineWidth', 1 );

   saveas(f,[datadir,filesep,'RhoB_depth_1.95_n=30_retraction.svg']); 

   save([datadir, filesep,'RhoB buildup_statistics_1.95_n=30.mat'],'FRET','FRET_arr','velocities', 'velocity_arr', 'myosin', 'myosin_arr'); %


