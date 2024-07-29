function correctBleachingExp2_cyto_ratio(position,datadir, fitpara_mRuby, fitpara_cyto)
% What it does
% Normalizes imRatio_raw for each frame by the median intensity per frame of the entire
% time series. 
% bleaching correction by division by linear bleaching function normalized
% by the median of its values. 

%remember to add in/remove 'position' tag for all the file names depending
%on if calling the stacked versions or not 

load([datadir,filesep,position, '_RatioData_raw.mat'],'imFRETOutline','im_mRuby_raw', 'im_cyto_raw');
load([datadir,filesep,position, '_Bleach_raw.mat']);
timepts=1:length(im_mRuby_raw);

corr_m=feval(fitpara_mRuby,timepts);
corr_norm_m=corr_m./median(corr_m);
normfact_m=nanmedian(bleach_raw_mRuby);

corr_cyto=feval(fitpara_cyto,timepts);
corr_norm_cyto=corr_m./median(corr_cyto);
normfact_cyto=nanmedian(bleach_raw_cyto)
 
 ezrin_ratio = cell(1,size(im_mRuby_raw,2));


for frameNum=1:length(im_mRuby_raw)
     
  
mRuby=im_mRuby_raw{frameNum}./(normfact_m*corr_norm_m(frameNum));
cyto = im_cyto_raw{frameNum}./(normfact_cyto*corr_norm_cyto(frameNum));    
 %bounds_mruby =[(prctile(mRuby{frameNum}, 5, 'all')),prctile(mRuby{frameNum}, 95, 'all')]
 %bounds_cyto =[(prctile(cyto{frameNum}, 5, 'all')),prctile(cyto{frameNum}, 95, 'all')]


  ezrin_ratio{frameNum} = mRuby ./ cyto; 
 ezrin_ratio{frameNum}(isnan(ezrin_ratio{frameNum}))=-100; 
  %bounds =[(prctile(ezrin_ratio{frameNum}, 5, 'all')),prctile(ezrin_ratio{frameNum}, 95, 'all')]

 
 colorRange2=[-2 4];
     
  
    tempmRubyforstack=ratio2RGB(ezrin_ratio{frameNum},colorRange2);%Cdc42'
    
    
    imwrite(tempmRubyforstack,[datadir,filesep,position, '_ezrin_[-2,4]','.tif'],'WriteMode','append','Compression','none');
   disp(num2str(frameNum)); 
%    
   
   
end

%load([datadir,filesep,position, '_RatioData_raw.mat'],'maskFinal','cellCoors');
save([datadir,filesep,position, '_CytoRatioData.mat'],'-v7.3','ezrin_ratio');
end    
