% immunostain compare to specific calculated ratio data 
%SM April 2024 
clc;
clear; 
ezrin_high_norm_array_95 = []; 
ezrin_high_norm_array_histo = []; 

pERM_high_norm_array_95 = []; 
pERM_high_norm_array_histo = []; 

%this must be run for both WT and T567 data! 
for i = 1:29 % number of cells you have in a single trial 
    
 
root = strcat('cropped data sites', num2str(i)); 
load([root, filesep,'output', filesep, 'ezrin_data', filesep, 'CytoRatioData.mat']); 
load([root, filesep,'output', filesep, 'RatioData_raw.mat']); 

datadir = ([root, filesep,'output', filesep, 'ezrin_data',filesep,'hotspot data']); 
if ~exist(datadir)
    mkdir(datadir)
end 
 pERM_raw=double(imread([root,filesep, 'pERM.tif']));
 mask = maskFinal{1,1}; 
 pERM_raw(~mask)=nan;
 pERM_high_95=pERM_raw; 
 pERM_high_histo=pERM_raw; 
 ezrin_high_95= ezrin_ratio{1}; 
 ezrin_high_histo= ezrin_ratio{1}; 
 
 %% %% section to segment based on pERM histogram 
 
 [f,xi]=ksdensity(pERM_raw(:)); 
% figure; plot(xi,f); hold on;  

[pks,locs]=findpeaks(f,xi);

log_pks=pks>0.0001; 
pks=pks(log_pks); locs=locs(log_pks);

   %     plot((locs),pks,'k^','markerfacecolor',[1 0 0]);

 
 
x_bgMax=locs(1,1); % picks the first peak (x-value of first peak)
[~,ind]=find(f>(0.5*pks(1)),1); % returns the first value of f greater than 0.5% of its max
x_1pct=xi(ind); % returns the corresponding intensity value
bgWidth=x_bgMax; % estimates the width of the background peak

  %    xline(x_bgMax+bgWidth); % to see how effective the estimation above of bg width is 



threshSeg=(x_bgMax+(8)*bgWidth);% 8 chosen as arbitrary threshold 


 
 %% continued section to segment hotspots with percentile thresholding option
 % 2 methods: histo uses > 8x bg intensity, while 95 uses > 95% percentile
 % values. You choose which you desire 
 
 bounds_pERM =[(prctile(pERM_raw, 5, 'all')),prctile(pERM_raw, 90, 'all')]; 
 
 pERM_high_95(pERM_raw <bounds_pERM(2)) = NaN; 
 ezrin_high_95(pERM_raw <bounds_pERM(2)) = NaN; 
 pERM_high_histo(pERM_raw <threshSeg) = NaN; 
 ezrin_high_histo(pERM_raw <threshSeg) = NaN; 
 
 
 
 stats = struct(); 
 stats.ezrin_mean = nanmean(ezrin_ratio{1,1}, 'all');   
 stats.ezrin_mean_high_95 = nanmean(ezrin_high_95, 'all'); 
 stats.ezrin_mean_high_histo = nanmean(ezrin_high_histo, 'all'); 
 stats.ezrin_high_norm_95 = stats.ezrin_mean_high_95/stats.ezrin_mean; 
  stats.ezrin_high_norm_histo = stats.ezrin_mean_high_histo/stats.ezrin_mean; 
 
 
stats.pERM_mean=nanmean(pERM_raw, 'all');  
stats.pERM_mean_high_95 = nanmean(pERM_high_95, 'all');  
stats.pERM_mean_high_histo = nanmean(pERM_high_histo, 'all');  
stats.pERM_high_norm_95 =  stats.pERM_mean_high_95/stats.pERM_mean; 
stats.pERM_high_norm_histo =  stats.pERM_mean_high_histo/stats.pERM_mean; 
 
  f1 = figure; 
 
  imagesc(pERM_high_95, [0 100]);
  hold off; 
   
  f2 = figure; 
  imagesc(pERM_high_histo, [0 100]);
  hold off; 
%   
  
 
  f3 = figure; 
 imagesc(ezrin_high_95, [0 2]);
 hold off; 
 
  f4= figure; 
 imagesc(ezrin_high_histo, [0 2]);
 hold off; 
 

 %% 
  saveas(f1,strcat(datadir, filesep,'pERM_high_90.fig'))
 saveas(f2,strcat(datadir, filesep,'pERM_high_histo.fig'))
  saveas(f3,strcat(datadir, filesep,'ezrin_high_90.fig'))
  saveas(f4,strcat(datadir, filesep,'ezrin_high_histo.fig'))


 save([datadir, filesep, 'pERM_ezrin_statistics.mat'], 'maskFinal', 'pERM_raw', 'pERM_high_95', 'pERM_high_histo' , 'ezrin_ratio', 'ezrin_high_histo','ezrin_high_95', 'stats', 'threshSeg');  
ezrin_high_norm_array_95 = [ezrin_high_norm_array_95, stats.ezrin_high_norm_95]; 
ezrin_high_norm_array_histo = [ezrin_high_norm_array_histo, stats.ezrin_high_norm_histo]; 
pERM_high_norm_array_95 = [pERM_high_norm_array_95, stats.pERM_high_norm_95]; 
pERM_high_norm_array_histo = [pERM_high_norm_array_histo, stats.pERM_high_norm_histo]; 
end 

 overalldata = struct(); 
 
overalldata.ezrin_average_90= mean(ezrin_high_norm_array_95);
overalldata.ezrin_average_histo = mean(ezrin_high_norm_array_histo);
overalldata.pERM_average_90 = mean(pERM_high_norm_array_95);
overalldata.pERM_average_histo = mean(pERM_high_norm_array_histo);

overalldata.ezrin_90= ezrin_high_norm_array_95; 
overalldata.ezrin_histo = ezrin_high_norm_array_histo; 
overalldata.pERM_90 = pERM_high_norm_array_95;
overalldata.pERM_histo = pERM_high_norm_array_histo;

 save(['data location',filesep,  'hotspot averages.mat'], 'overalldata');  


%% for comparing WT vs T567
clear; 
load('t567 data saved in above section location\hotspot averages.mat'); 

t567overall = overalldata; 

load('WT data saved in above section location\hotspot averages.mat'); 

WToverall = overalldata; 

pERM = struct; 
pERM.WT = WToverall.pERM_histo; 
pERM.t567 = t567overall.pERM_histo; 

violinplot(pERM); 
hold off; 

f2 = figure; 
ylabel('Normalized ezrin intensity'); 
title('pERM hotspots: 8x histo'); 
ylim([0.8 1.8]); 
ezrin = struct; 
ezrin.WT = WToverall.ezrin_histo; 
ezrin.t567 = t567overall.ezrin_histo; 

violinplot(ezrin); 


%% combine trials together: 

clear; 
load('t567 T1 filepath\cropped\hotspot averages.mat'); 

t567ezrin = overalldata.ezrin_histo; 
t567pERM = overalldata.pERM_histo; 

load('t567 T2 filepath\cropped\hotspot averages.mat'); 

t567ezrin = [t567ezrin, overalldata.ezrin_histo]; 
t567pERM = [t567pERM, overalldata.pERM_histo];  
 
load('WT T1 filepath\cropped\hotspot averages.mat'); 

WTezrin = overalldata.ezrin_histo; 
WTpERM = overalldata.pERM_histo; 

load('WT T2 filepath\cropped\hotspot averages.mat'); 

WTezrin = [WTezrin, overalldata.ezrin_histo]; 
WTpERM = [WTpERM, overalldata.pERM_histo];  

s1 = struct; 
s1.WT = WTezrin; 
s1.t567 = t567ezrin; 
ylabel('Norm. ezrin ratio intensity'); 
title('Ezrin ratio in thresholded hotspots'); 
violinplot(s1); 
ylim([0.6 2.4]); 
axis square; 

ranksum(t567ezrin, WTezrin) 



