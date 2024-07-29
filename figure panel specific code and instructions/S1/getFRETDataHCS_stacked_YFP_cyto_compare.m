function getFRETDataHCS_stacked_YFP_cyto_compare(cellNum,rawdir,datadir)
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
if ~exist([datadir,filesep,'RatioData_raw.mat'])
%%%%%% Set up

load([rawdir,filesep,'alignment parameters pX pY.mat'],'pX','pY');

%%%%%% Call background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
CFPbg_raw=double(imread([rawdir,filesep,'AVG_bgCFP.tif']));
FRETbg_raw=double(imread([rawdir,filesep,'AVG_bgFRET.tif']));
cytobg_raw=double(imread([rawdir,filesep,'AVG_bgmRuby.tif']));
bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);

 [m,n]=size(cytobg_raw);
 bg3(:,:,1)=cytobg_raw;bg3(:,:,2)=zeros(m,n);
 bg4=dualviewAlignFromFittedSurface(bg3,pX,pY,binning);
 cytobg=bg4(:,:,1);

cellPath=strcat('cell_',num2str(cellNum));
cellFiles=getFilenames([rawdir],'.tif');

 CFP_stack=double(readTIFFstack([rawdir,filesep,cellFiles{4}]));
 FRET_stack=double(readTIFFstack([rawdir,filesep,cellFiles{5}]));
 cyto_stack=double(readTIFFstack([rawdir,filesep,cellFiles{6}]));
 
  
%%%%%% Loop through frames
imRatio_raw={};maskFinal={};cellCoors={};   im_cyto_raw={};
for frameNum=1:size(CFP_stack,3)
    disp(num2str(frameNum));
    imCFP_raw=CFP_stack(:,:,frameNum);
    imFRET_raw=FRET_stack(:,:,frameNum);
    imcyto_raw=cyto_stack(:,:,frameNum);

    %%%%%% Align CFP/FRET images
    imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
    imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
    imCFP_raw=imaligned(:,:,1);
    imFRET_raw=imaligned(:,:,2);
    
     imstack2(:,:,1)=imcyto_raw; imstack2(:,:,2)=zeros(m,n);
     imaligned2=dualviewAlignFromFittedSurface(imstack2,pX,pY,1);
     imcyto_raw=imaligned2(:,:,1);
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imCFP_raw+imFRET_raw);
    bgmask_cyto=getBGMask(imcyto_raw); % -- having issues here with

    imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
    imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
    imcytobg=subBG(imcyto_raw,bgmask_cyto,cytobg);
    %%%%%% Get mask from raw FRET image
    s = cellNum; 
   [mask cellCoorsTemp]=getCellMaskCyto_2_stacked((imFRETbg),1500,frameNum,s); % see here if its better taking away the imCFPbg
   
    maskFinal{frameNum}=mask;
   cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Detrmine ratio
    imFRETbg(~mask)=nan;
    imFRET=ndnanfilter(imFRETbg,fspecial('disk',3),'replicate');
    imCFPbg(~mask)=nan;
    imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');
    imRatioTemp=imFRET./imCFP;
    imRatioTemp(~mask)=nan;
    imRatio_raw{frameNum}=imRatioTemp;
    imYFP_raw{frameNum}=imFRETbg; 
    
     imcytobg(~mask)=nan;
    im_cyto_raw{frameNum}=imcytobg;

   
    %%%%%% Generate and write files for raw ratio and outlined objects - optional
    %%%%%% 
   
    imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,mask);

   imwrite(imFRETOutline{frameNum},[datadir,filesep,'Outline_prelim.tif'],'WriteMode','append','Compression','none');
   
end
%%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay

for frameNum= 1:length(imRatio_raw)
    bleach_YFP(frameNum)=nanmean(vect(imYFP_raw{frameNum})); 
    bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
    bleach_raw_cyto(frameNum)= nanmean(vect(im_cyto_raw{frameNum}));
end
save([datadir,filesep,'RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw','imFRETOutline','imYFP_raw','im_cyto_raw, '-v7.3');
save([datadir,filesep,'Bleach_raw.mat'],'bleach_raw' ,  'bleach_YFP', 'bleach_raw_cyto'); 
end

