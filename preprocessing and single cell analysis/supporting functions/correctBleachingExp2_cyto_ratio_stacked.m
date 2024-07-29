function correctBleachingExp2_cyto_ratio_stacked(datadir, fitpara_mRuby, fitpara_cyto)
% What it does
% Normalizes imRatio_raw for each frame by the median intensity per frame of the entire
% time series. 
% bleaching correction by division by linear bleaching function normalized
% by the median of its values. 


load([datadir,filesep, 'RatioData_raw.mat'],'imFRETOutline','im_mRuby_raw', 'im_cyto_raw');
load([datadir,filesep,'Bleach_raw.mat']);
timepts=1:length(im_mRuby_raw); 

corr_m=feval(fitpara_mRuby,timepts);
corr_norm_m=corr_m./median(corr_m);
normfact_m=nanmedian(bleach_raw_mRuby);

corr_cyto=feval(fitpara_cyto,timepts);
corr_norm_cyto=corr_m./median(corr_cyto);
normfact_cyto=nanmedian(bleach_raw_cyto);

 ezrin_ratio = cell(1,size(im_mRuby_raw,2)); 


for frameNum=1:length(im_mRuby_raw) 
     
  
mRuby{frameNum}=im_mRuby_raw{frameNum}./(normfact_m*corr_norm_m(frameNum));
cyto{frameNum} = im_cyto_raw{frameNum}./(normfact_cyto*corr_norm_cyto(frameNum));    



 ezrin_ratio{frameNum} = mRuby{frameNum} ./ cyto{frameNum}; 
 
 ezrin_ratio_new{frameNum}(isnan(ezrin_ratio{frameNum}))=0; 
  bounds =[(prctile(ezrin_ratio{frameNum}, 1, 'all')),prctile(ezrin_ratio{frameNum}, 99, 'all')]; 

 
   colorRange=[0 3]; 
    tempmRubyforstack=ratio2RGB(mRuby{frameNum},colorRange);
    tempmRubyforstack(tempmRubyforstack < 0) = 0; 
    imwrite(tempmRubyforstack,[datadir,filesep, 'mRuby_[0,3]','.tif'],'WriteMode','append','Compression','none');
 
    
    colorRange2=[0.5 1.5];
    % here essentially removing extreme outlier values by setting them to
    % te 1st and 99th precentile 
   ezrin_ratio{frameNum}(ezrin_ratio{frameNum}< bounds(1,1))=bounds(1,1);
   ezrin_ratio{frameNum}(ezrin_ratio{frameNum}> bounds(1,2))=bounds(1,2);

    tempmRubyforstack2=ratio2RGB(ezrin_ratio{frameNum},colorRange2);%Cdc42'
   imwrite(tempmRubyforstack2,[datadir,filesep, 'ezrin_[0.5 1.5]','.tif'],'WriteMode','append','Compression','none');
   disp(num2str(frameNum)); 
%    
   
   
end

%load([datadir,filesep,position, '_RatioData_raw.mat'],'maskFinal','cellCoors');
save([datadir,filesep, 'CytoRatioData.mat'],'-v7.3','ezrin_ratio','mRuby','cyto');
end    
