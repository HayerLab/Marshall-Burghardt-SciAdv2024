function correctBleachingExp2(position,fitpara,datadir, fitpara_mRuby)

% Normalizes imRatio_raw for each frame by the median intensity per frame of the entire
% time series. 
% bleaching correction by division by linear bleaching function normalized
% by the median of its values. 



load([datadir,filesep,position, '_RatioData_raw.mat'],'imRatio_raw','imFRETOutline','im_mRuby_raw');
load([datadir,filesep,position, '_Bleach_raw.mat']);
timepts=1:length(imRatio_raw);

 corr=feval(fitpara,timepts); 
 corr_norm=corr./median(corr); 
 normfact=nanmedian(bleach_raw); 


corr_m=feval(fitpara_mRuby,timepts);
corr_norm_m=corr_m./median(corr_m);
normfact_m=nanmedian(bleach_raw_mRuby);



 imRatio= cell(1,size(imRatio_raw,2));
 im_mRuby= cell(1,size(imRatio_raw,2));
 %RATIO_cyto = cell(1,size(imRatio_raw,2));
 %mRuby_cyto = cell(1,size(imRatio_raw,2));
for frameNum=1:length(imRatio_raw)
   
    imRatio{frameNum}=imRatio_raw{frameNum}./(normfact*corr_norm(frameNum));
 colorRange = [0.7 1.3]; 
 
     tempRATIOforstack=ratio2RGB( imRatio{frameNum},colorRange);
    
   
  imwrite(tempRATIOforstack,[datadir,filesep,position,  '_Rho-FRET','.tif'],'WriteMode','append','Compression','none');
 
% imwrite(imFRETOutline{frameNum},[datadir,filesep,'_Outline','.tif'],'WriteMode','append','Compression','none');

 

 
end

for frameNum=1:length(im_mRuby_raw)
     
  
 im_mRuby{frameNum}=im_mRuby_raw{frameNum}./(normfact_m*corr_norm_m(frameNum));
   
     colorRange2=[0 3];
     
    tempmRubyforstack=ratio2RGB(im_mRuby{frameNum},colorRange2);
   
    tempmRubyforstack(tempmRubyforstack < 0) = 0; 
    
    imwrite(tempmRubyforstack,[datadir,filesep,position, '_mRuby_[0,3]','.tif'],'WriteMode','append','Compression','none');
   disp(num2str(frameNum)); 
   
   
   
end

load([datadir,filesep,position, '_RatioData_raw.mat'],'maskFinal','cellCoors');
save([datadir,filesep,position, '_RatioData.mat'],'-v7.3','maskFinal','cellCoors','imRatio','imFRETOutline','im_mRuby');
end    
