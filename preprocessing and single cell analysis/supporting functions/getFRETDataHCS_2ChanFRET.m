function getFRETDataHCS_2ChanFRET(position,bgdir,rawdir,datadir,threshold)
%function getFRETDataHCS( row,col,site )
% image processing for FRET data analysis 
% row=3;col=3;site=2;
% Input: CFP/FRET images, background images,  and alignment parameters for
% CFP/FRET images. 
%
% Generates first a mask for background subtraction, subtracts background,
% and calculates a refined mask before taking a ratio between background-
% subtracted CFP and FRET images. Bleaching correction assumes constant 
% mean FRET values per frame across the entire time-lapse series.
%
% Output: 
%
%   maskFinal   Cell array with binary mask for each frame
%
%   cellCoors   Cell array with centroid x,y coordinates of detected
%               objects and object size in pixels
%
%   imRatioraw  Cell array with masked raw ratio values, before bleaching
%               correction
%
%   imRatio     Cell array with bleaching-corrected masked ratio values
%
%   imFRETOUtline_row_col.tif 
%               RGB tif stack of FRET channel images with outlined detected
%               objects
%
%   imRatio_row_col_colorRangeLow_colorRangeHigh.tif
%               RGB tif stack with false-colored (jet-black) scaled ratio
%               images, bleaching corrected.
%
% Used non built-in subfunctions: 
%
%   dualviewAlignFromFittedSurface
%   getBGMask
%   subBG
%   getCellMask
%   DrawMaskOutline
%   ratio2RGB
%   vect
%
% Arnold, 24 May 2015
key=[datadir,filesep,position,'_RatioData_raw.mat'];
if ~exist(key)
%%%%%% Set up

load([bgdir,filesep,'alignment parameters pX pY_day2.mat'],'pX','pY');

%%%%%% Call background images
 binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
CFPbg_raw=double(imread([bgdir,filesep,'AVG_bgCFP.tif']));
FRETbg_raw=double(imread([bgdir,filesep,'AVG_bgFRET.tif']));
 
bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
bg2=dualviewAlignFromFittedSurfaceCrop(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);

 
CFP_files=dir([rawdir,filesep,position,'_CFP_*']);
%%%%%% Loop through frames
imRatio_raw={};maskFinal={};cellCoors={}; imFRETOutline = {};   % im_mRuby_raw={};
for frameNum=1:length(CFP_files)
   
  




    disp([position,'__',num2str(frameNum)]);
    imCFP_raw=double(imread([rawdir,filesep,position,'_CFP_',num2str(frameNum),'.tif']));
    imFRET_raw=double(imread([rawdir,filesep,position,'_FRET_',num2str(frameNum),'.tif']));
  
   

    %%%%%% Align CFP/FRET images
    imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
    imaligned=dualviewAlignFromFittedSurfaceCrop(imstack,pX,pY,1);
    imCFP_raw=imaligned(:,:,1);
    imFRET_raw=imaligned(:,:,2);
%     
    
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imCFP_raw+imFRET_raw);
    imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
    imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
 
   %%%%%% Get mask from raw FRET image
    [mask cellCoorsTemp]=getCellMaskCyto_2(imFRETbg,1500, frameNum, position, threshold); %+imCFPbg
 
    maskFinal{frameNum}=mask;
    cellCoors{frameNum}=cellCoorsTemp;
   
    %%%%%% Detrmine ratio
    imFRETbg(~mask)=nan;
    imFRET=ndnanfilter(imFRETbg,fspecial('disk',3),'replicate');
    imCFPbg(~mask)=nan;
    imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');
    imRatio_raw{frameNum}=imFRET./imCFP;
    imRatio_raw{frameNum}(~mask)=nan;
    
    
    %%%%%% Generate and write files for raw ratio and outlined objects
   % tempRATIO=ratio2RGB(imRatio_raw{frameNum},colorRange);
    imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,mask);
    
    imwrite(imFRETOutline{frameNum},[datadir,filesep,position,'_Outline_prelim.tif'],'WriteMode','append','Compression','none');
end
%%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay
bleach_1=nanmean(vect(imRatio_raw{1}));
bleach_raw = {};  
for frameNum=1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
 
end
save([datadir,filesep,position,'_RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw','imFRETOutline','-v7.3');  

save([datadir,filesep,position,'_Bleach_raw.mat'],'bleach_raw'); % 'CFP_raw','FRET_raw'); 
end  

