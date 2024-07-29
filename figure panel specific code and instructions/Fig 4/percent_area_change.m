%% percent area change 
% script to calculate changes in cell size pre/post drug treatment 3
clc; clear; 

 root1 = 'F:\example dataset'; 
 
 %% 
 
 
 change_control = []; 
 change_Y2 = []; 
 change_ML7 = []; 
 change_Y2_ML7=[]; 
 
 k=0;

 %control sites
for row=1:2    
    for col=1
        for site=1:4


            k=k+1;
            position_cntrl{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end

  k=0;
%drug condition 1 sites 
for row=1:2   
    for col=1
        for site=5:8


            k=k+1;
            position_Y2{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end

 k=0;
 %drug condition 2 sites 
for row=1:2 
    for col=1
        for site=9:12

            k=k+1;
            position_ML7{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end
 
 k=0;
 %drug position 3 sites 
for row=1:2
    for col=1
        for site=13:16

            k=k+1;
            position_Y2_ML7{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end
 
%% control
percent_change_cntrl = []; 

for i = 1: length(position_cntrl)
    load([root1, filesep,strcat(position_cntrl{i}, '_RatioData_raw.mat')]); 
    mask_pre = []; 
    mask_post = []; 
   mask_pre = maskFinal{1,5}(maskFinal{1,5} ==1); 
   mask_post = maskFinal{1,35}(maskFinal{1,35} ==1); 
   
   percent_change = ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; 
   
   percent_change_cntrl = [percent_change_cntrl,percent_change]; 

end 
 
%% drug condition 1 
percent_change_Y2 = []; 


for i = 1: length(position_Y2)
    mask_pre =[]; 
    mask_post =[]; 
    load([root1, filesep,strcat(position_Y2{i}, '_RatioData_raw.mat')]); 
    
   mask_pre = maskFinal{1,5}(maskFinal{1,5} ==1); 
   mask_post = maskFinal{1,35}(maskFinal{1,35} ==1); 
   
   percent_change = ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; 
   
   percent_change_Y2 = [percent_change_Y2,percent_change]; 

end 

%% drug condition 2 
percent_change_ML7 = []; 


for i = 1: length(position_ML7)
    mask_pre =[]; 
    mask_post =[]; 
    load([root1, filesep,strcat(position_ML7{i}, '_RatioData_raw.mat')]); 
    
   mask_pre = maskFinal{1,5}(maskFinal{1,5} ==1); 
   mask_post = maskFinal{1,35}(maskFinal{1,35} ==1); 
   
   percent_change = ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; 
   
   percent_change_ML7 = [percent_change_ML7,percent_change]; 

end 
 
%% drug condition 3 
percent_change_Y2_ML7 = []; 


for i = 1: length(position_Y2_ML7)
    mask_pre =[]; 
    mask_post =[]; 
    load([root1, filesep,strcat(position_Y2_ML7{i}, '_RatioData_raw.mat')]); 
    
   mask_pre = maskFinal{1,5}(maskFinal{1,5} ==1); 
   mask_post = maskFinal{1,35}(maskFinal{1,35} ==1); 
   
   percent_change = ((size(mask_post,1)-size(mask_pre,1))/size(mask_pre,1))*100; 
   
   percent_change_Y2_ML7 = [percent_change_Y2_ML7,percent_change]; 

end 

%% 
 save('save where you want', 'percent_change_cntrl','percent_change_Y2', 'percent_change_ML7', 'percent_change_Y2_ML7');  
 
 
 %% load multiple trials of the above code to plot 
 clear; 
 cntrl =[];
 Y2 =[]; 
 ML7 =[]; 
 Y2_ML7 = []; 
 
 load('trial 1 path.mat'); 
 
 cntrl = [cntrl, percent_change_cntrl]; 
 Y2 =[Y2, percent_change_Y2]; 
 ML7 =[ML7, percent_change_ML7]; 
 Y2_ML7 = [Y2_ML7, percent_change_Y2_ML7];
 
 load('trial 2 path.mat'); 
 
 cntrl = [cntrl, percent_change_cntrl]; 
 Y2 =[Y2, percent_change_Y2]; 
 ML7 =[ML7, percent_change_ML7]; 
 Y2_ML7 = [Y2_ML7, percent_change_Y2_ML7];
 
  
 load('trial 3 path.mat'); 
 
 cntrl = [cntrl, percent_change_cntrl]; 
 Y2 =[Y2, percent_change_Y2]; 
 ML7 =[ML7, percent_change_ML7]; 
 Y2_ML7 = [Y2_ML7, percent_change_Y2_ML7];
 
 %% 
 s = struct; 



s.Y2_ML7 = Y2_ML7; 
s.ML7 = ML7; 
s.DY2 = Y2; 
s.CNTRL = cntrl; 

f= figure; 

grid on; 
violinplot(s); 

yticks(-100:10:40); 
ylim([-100, 40]); 
daspect([1 1 1]); 

% 
  y = [ cntrl, Y2, ML7, Y2_ML7]; 
  group = repelem(1:4, 1, [numel(cntrl),numel(Y2),numel(ML7), numel(Y2_ML7)]);
 
  %% riun statistics 
 [p, anova, stats] = anova1(y, group); 
 
 [compare]= multcompare(stats) 

 save('fielpath\compiled.mat', 'compare', 'Y2', 'ML7', 'cntrl', 'Y2_ML7'); 
 
 
 
 
 
 