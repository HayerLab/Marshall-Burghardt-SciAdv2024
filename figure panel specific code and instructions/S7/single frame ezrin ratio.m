% for single frame analysis of ezrin ratio

clc; clear; 
for j=1:20 % number of cells you wish to analyze 

close all; 
root = ['file location']; 
datadir=[root, filesep, 'ezrin data']; 
if ~exist(datadir)
    mkdir(datadir)
end 

load([root, filesep, 'RatioData_raw.mat'],'im_mRuby_raw', 'im_cyto_raw'); 
load([root, filesep, 'Bleach_raw.mat']); 


timepts=1:length(im_mRuby_raw);

frames= 1; 
mRuby_predrug = []; 
cyto_predrug = []; 

for i = 1:frames
  
    mRuby_predrug = [mRuby_predrug, nanmean(im_mRuby_raw{1, i}, 'all')]; 
    cyto_predrug = [cyto_predrug, nanmean(im_cyto_raw{1, i}, 'all')]; 
   
end 
 
mRuby_mean  = sum(mRuby_predrug)/frames; 
cyto_mean = sum(cyto_predrug)/frames;


ezrin_ratio = cell(1,size(im_mRuby_raw,2));
 ezrin_mean = []; 


for frameNum=1:length(im_mRuby_raw)
     
  
mRuby{frameNum}=im_mRuby_raw{frameNum}   ./mRuby_predrug(1, frameNum); 
cyto{frameNum} = im_cyto_raw{frameNum}   ./cyto_predrug(1, frameNum);  

 ezrin_ratio{frameNum} = mRuby{frameNum} ./ cyto{frameNum}; 

 bounds =[(prctile(ezrin_ratio{frameNum}, 1, 'all')),prctile(ezrin_ratio{frameNum}, 99, 'all')];
   ezrin_ratio{frameNum}(ezrin_ratio{frameNum}< bounds(1,1))=bounds(1,1);
    ezrin_ratio{frameNum}(ezrin_ratio{frameNum}> bounds(1,2))=bounds(1,2);
 ezrin_mean(1,frameNum) = nanmean(ezrin_ratio{1, frameNum}, 'all'); 

 ezrin_ratio_forRGB{frameNum}=ezrin_ratio{frameNum};
 ezrin_ratio_forRGB{frameNum}(isnan(ezrin_ratio_forRGB{frameNum}))=-100; 
 
 ezrin_ratio_norm{frameNum} = ezrin_ratio{frameNum};   % ./ezrin_mean(1, frameNum); 

   colorRange=[0 2]; 
    tempmRubyforstack=ratio2RGB(mRuby{frameNum},colorRange);
    tempmRubyforstack(tempmRubyforstack < 0) = 0; 
    imwrite(tempmRubyforstack,[datadir,filesep, 'mRuby_[0,2]','.tif'],'WriteMode','append','Compression','none');
 
    
    colorRange2=[0 2];
    colorRange3=[-5 6];
    

    tempmRubyforstack2=ratio2RGB(ezrin_ratio_forRGB{frameNum},colorRange2);%Cdc42'
   imwrite(tempmRubyforstack2,[datadir,filesep, 'ezrin_[0 2]','.tif'],'WriteMode','append','Compression','none');
   


   
end

save([datadir,filesep, 'CytoRatioData.mat'],'-v7.3','ezrin_ratio', 'ezrin_mean', 'ezrin_ratio_norm'); 
disp(j); 
clear; 
end 
