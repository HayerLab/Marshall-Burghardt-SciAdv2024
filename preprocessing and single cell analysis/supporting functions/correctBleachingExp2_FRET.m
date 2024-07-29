function correctBleachingExp2_FRET(position,fitpara,datadir)
% What it does
% Normalizes imRatio_raw for each frame by the median intensity per frame of the entire
% time series. 
% bleaching correction by division by linear bleaching function normalized
% by the median of its values. 

%remember to add in/remove 'position' tag for all the file names depending
%on if calling the stacked versions or not 

load([datadir,filesep,position, '_RatioData_raw.mat'],'imRatio_raw');
load([datadir,filesep,position, '_Bleach_raw.mat']);
timepts=1:length(imRatio_raw);

 corr=feval(fitpara,timepts); 
 corr_norm=corr./median(corr); 
 normfact=nanmedian(bleach_raw);
 



 imRatio= cell(1,size(imRatio_raw,2));
 
for frameNum=1:length(imRatio_raw)
   
    imRatio{frameNum}=imRatio_raw{frameNum}./(normfact*corr_norm(frameNum));
 
colorRange = [0.7 1.3]; 
tempRATIOforstack=ratio2RGB( imRatio{frameNum},colorRange);
    
   
  imwrite(tempRATIOforstack,[datadir,filesep,position,  '_Rho-FRET','.tif'],'WriteMode','append','Compression','none');


end



load([datadir,filesep,position, '_RatioData_raw.mat'],'maskFinal','cellCoors','imFRETOutline');
save([datadir,filesep,position, '_RatioData.mat'],'-v7.3','maskFinal','cellCoors','imRatio','imFRETOutline');
end    
