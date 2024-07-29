function correctBleachingExp2_stacked_YFP_cyto(datadir)
% What it does
% Normalizes imRatio_raw for each frame by the median intensity per frame of the entire
% time series. 
% bleaching correction by division by linear bleaching function normalized
% by the median of its values. 


load([datadir,filesep, 'RatioData_raw.mat'],'imRatio_raw','imFRETOutline','imYFP_raw', 'im_cyto_raw');
load([datadir,filesep,'Bleach_raw.mat']);




 bleach_YFP_all=bleach_YFP/nanmedian(bleach_YFP);

 

timepts=1:length(imRatio_raw);

im_YFP = cell(1,size(imYFP_raw,2));


for frameNum=1:length(imYFP_raw)
     
  
 im_YFP{frameNum}=imYFP_raw{frameNum}./(bleach_YFP_all(1,frameNum)* bleach_YFP(1, frameNum));
 
    

    colorRange2=[0 2.5];
    tempforstack=ratio2RGB(im_YFP{frameNum},colorRange2);%Cdc42'
    tempforstack(tempforstack < 0) = 0.001; 
    
    imwrite(tempforstack,[datadir,filesep, 'YFP_[0,2.5]','.tif'],'WriteMode','append','Compression','none');
   disp(num2str(frameNum)); 
  
 
end

load([datadir,filesep, 'RatioData_raw.mat'],'maskFinal','cellCoors');
save([datadir,filesep, 'RatioData.mat'],'-v7.3','maskFinal','cellCoors','imRatio','imFRETOutline','im_YFP');
end    