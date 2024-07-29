%designed to measure extent of FRET activity preceding specific
%protrusion/retraction event 

%% data loading, initialization 



clc; clear; 
close all; 

root  = 'F:\example path'; 
event = 1; 
rawdir = ([root, filesep,'cropped',filesep, '1']); 
datadir =([root, filesep,'FRET buildup - retractions', filesep, num2str(event)]); 

if  ~exist(datadir)
    mkdir(datadir)
end 

depths =[3,6,10,15,20,25]; 

FRET_depths = cell(2,6);
%myosin_depths = cell(2,6); 

load([rawdir, filesep,'edge_vels', filesep,  'edge vel mapping_3',filesep,'Protrusion and FRET values.mat'],'protvalsWindowF'); 

for row = 1: size(protvalsWindowF,1)
        for col = 1:size(protvalsWindowF,2)-1
          temporary(row,col)=(protvalsWindowF(row,col)+protvalsWindowF(row,col+1))/2;
        end
end

protval_map = temporary; 

for i=1:size(depths, 2)
   
    load([rawdir, filesep,strcat( 'edge_vels', filesep, 'edge vel mapping_',num2str(depths(1,i))),filesep,'Protrusion and FRET values.mat'],'fretvalsF','myosinF')
    
    depth =  num2str(depths(1,i)); 
    
   FRET_depths{1,i}=depth;  
   FRET_depths{2,i}= fretvalsF(:,2:end-1);
   
%    myosin_depths{1,i}=depth;  
%    myosin_depths{2,i}= myosinF(:,2:end-1);
    
end 


load('CMAP_blue_grey_yellow.mat');
w=figure;
imagesc(protval_map,[-13,13]);title('Edge Velocity');
colormap(w,cmap);


%% specify rectangle dimensions and plot Rho activity vs myosin localization 
close all; 
f = figure; 
 
imagesc(protval_map,[-13,13]);title('Edge Velocity');
colormap(f,cmap);

hold on; 
%manually input these to specify rough ROI for protrusion/retraction event
%of interest 
coor_start = 140;
time_start=75;

coor_end = coor_start +15;
time_end = time_start+60; 

if time_end >= size(protval_map,2)
    time_end =size(protval_map,2);
    
end 

 
rectangle('Position',[time_start,coor_start,time_end-time_start,coor_end-coor_start],'LineWidth',2); 

hold off; 

f2 = figure; 


region = protval_map(coor_start:coor_end, time_start:time_end); 



% calculate average velocity in specified area 
for t = 1:size(region,2)
vel_avg(1,t) = mean(region(:,t)); 
end 


hold on; 
%specify how many frames pre and post t = 0 (optional) 
x = [-15:1:30]; 



%find timept of max acceleration (pos or neg depending on analysis)
acceleration=zeros(1,size(vel_avg,2)-1);
for i=1:size(vel_avg,2)-1
    acceleration(1,i)=[vel_avg(1,i+1)-vel_avg(1,i)]; 
end 

% change this depending on if a protrusion or retraction 
for j = 1:size(vel_avg,2)
    %protrusion to retraction 
  zero_pt = (acceleration(1,j) <0  && vel_avg(1,j) >= 0 && vel_avg(1,j+1) <= 0); 
   %retraction to protrusion
   % zero_pt = (acceleration(1,j) >0  && vel_avg(1,j) <= 0 && vel_avg(1,j+1) >= 0); 
   if (zero_pt ==1 && j > 3)
       
       break;   
   end 

end 

if (abs(vel_avg(1,j+1)) < abs(vel_avg(1,j)))

    align_pt =j+1; 
    
else 
    align_pt = j; 

end 
%xline(0); 

%here ensure that these dimensions match the ones you set on line 97
vel_arr = vel_avg(1,align_pt-15:align_pt+30);
yyaxis left;
ylabel('Myosin'); 
xlabel('Timepoints');
%ylim([-15 15 ]); 

myosin_temp=cell(2,6); 
for z =1:6
myosin_temp{1,z}=depths(1,z); 
myosin_temp{2,z} = myosin_depths{2,z}(coor_start:coor_end,time_start+align_pt-1-15:time_start+align_pt-1+30); 
myosin_temp{2,z}(isnan(myosin_temp{2,z}))=0; 
for w=1:size(myosin_temp{2,z},2)
myosin_avg(1,w) =mean(myosin_temp{2,z}(:,w));     
end 

plot(x,myosin_avg); 
 end 
 
 
xline(0, '--'); 


yyaxis right;
ylabel('RhoA');
%ylim([0.5 1.5]); 

%yline(1, '--k'); 

FRET_temp=cell(2,6); 
for z =1:6
FRET_temp{1,z}=depths(1,z); 
FRET_temp{2,z} = FRET_depths{2,z}(coor_start:coor_end,time_start+align_pt-1-15:time_start+align_pt-1+30); 
FRET_temp{2,z}(isnan(FRET_temp{2,z}))=0; 
for w=1:size(FRET_temp{2,z},2)
FRET_avg(1,w) =mean(FRET_temp{2,z}(:,w));     
end 

plot(x,FRET_avg); 
 
end 



hold off; 
 %legend('','','0.975 um','1.95 um', '3.25 um', '4.88 um', '6.50 um','8.13 um','Location','southoutside'); 



 saveas(f, [datadir,filesep,'Vel Map.png']); 
 saveas(f2,[datadir,filesep,'retraction']); 
% 
 save(strcat(datadir,'\','retraction_statistics_Rho vs_myosin.mat'),'time_start','time_end','coor_start','coor_end','FRET_temp','vel_arr','myosin_temp', 'align_pt');
 %
%% example of plotting Rho vs edge velocity instead

vel_arr = vel_avg(1,align_pt-15:align_pt+30); 

f3 = figure; 

x =-15:30; 



yyaxis right;
ylabel('RhoA');
%ylim([0.8 1.5]); 
FRET_temp=cell(2,6); 
for z =1:6
FRET_temp{1,z}=depths(1,z); 
FRET_temp{2,z} = FRET_depths{2,z}(coor_start:coor_end,time_start+align_pt-1-10:time_start+align_pt-1+30); 
FRET_temp{2,z}(isnan(FRET_temp{2,z}))=0; 
for w=1:size(FRET_temp{2,z},2)
FRET_avg(1,w) =mean(FRET_temp{2,z}(:,w));     
end 

plot(x,FRET_avg); 
 
end 

xline(0, '--'); 

yyaxis left; 
ylabel('Edge Velocity')
plot(x, vel_arr); 
 
  saveas(f3,[datadir,filesep,'retraction']); 
% % 
  save(strcat(datadir,'\','retraction_statistics.mat'),'time_start','time_end','coor_start','coor_end','FRET_temp','vel_arr','myosin_temp', 'align_pt'); 
% 
