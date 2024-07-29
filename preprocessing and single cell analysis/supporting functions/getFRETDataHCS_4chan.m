function getFRETDataHCS_4chan(position,bgdir,rawdir,datadir,threshold)
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

load([bgdir,filesep,'alignment parameters pX pY.mat'],'pX','pY');

%%%%%% Call background images
 binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
CFPbg_raw=double(imread([bgdir,filesep,'AVG_bgCFP.tif']));
FRETbg_raw=double(imread([bgdir,filesep,'AVG_bgFRET.tif']));
 mRubybg_raw=double(imread([bgdir,filesep,'AVG_bgmRuby.tif']));
cytobg_raw = double(imread([bgdir, filesep, 'AVG_bgCyto.tif'])); 

bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
bg2=dualviewAlignFromFittedSurfaceCrop(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);

 [m,n]=size(mRubybg_raw);
 bg3(:,:,1)=mRubybg_raw;bg3(:,:,2)=zeros(m,n);
 bg4=dualviewAlignFromFittedSurfaceCrop(bg3,pX,pY,binning);
 mRubybg=bg4(:,:,1);
 
  [j,k]=size(cytobg_raw);
 bg5(:,:,1)=cytobg_raw;bg5(:,:,2)=zeros(j,k);
 bg6=dualviewAlignFromFittedSurfaceCrop(bg5,pX,pY,binning);
 cytobg=bg6(:,:,1);

CFP_files=dir([rawdir,filesep,position,'_CFP_*']);
%%%%%% Loop through frames
imRatio_raw={};maskFinal={};cellCoors={}; imFRETOutline = {};  im_mRuby_raw={}; im_cyto_raw = {}; 
for frameNum=1:length(CFP_files)
   
  




    disp([position,'__',num2str(frameNum)]);
    imCFP_raw=double(imread([rawdir,filesep,position,'_CFP_',num2str(frameNum),'.tif']));
    imFRET_raw=double(imread([rawdir,filesep,position,'_FRET_',num2str(frameNum),'.tif']));
    imRuby_raw=double(imread([rawdir,filesep,position,'_mRuby_',num2str(frameNum),'.tif']));
    imcyto_raw=double(imread([rawdir,filesep,position,'_cyto_',num2str(frameNum),'.tif']));
   
   

    %%%%%% Align CFP/FRET images
    imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
    imaligned=dualviewAlignFromFittedSurfaceCrop(imstack,pX,pY,1);
    imCFP_raw=imaligned(:,:,1);
    imFRET_raw=imaligned(:,:,2);
%     
     imstack2(:,:,1)=imRuby_raw; imstack2(:,:,2)=zeros(m,n);
     imaligned2=dualviewAlignFromFittedSurfaceCrop(imstack2,pX,pY,1);
    imRuby_raw=imaligned2(:,:,1);
    
    imstack3(:,:,1)=imcyto_raw; imstack3(:,:,2)=zeros(j,k);
     imaligned3=dualviewAlignFromFittedSurfaceCrop(imstack3,pX,pY,1);
    imcyto_raw=imaligned3(:,:,1);
    
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imCFP_raw+imFRET_raw);
    bgmask_mRuby=getBGMask(imRuby_raw);
    bgmask_cyto=getBGMask(imcyto_raw);
    imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
    imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
    imRubybg=subBG(imRuby_raw,bgmask_mRuby,mRubybg);
     imcytobg=subBG(imcyto_raw,bgmask_cyto,cytobg);
    %%%%%% Get mask from raw FRET image
    [mask cellCoorsTemp]=getCellMaskCyto_2(imFRETbg,5000, frameNum, position, threshold); 
   
    

    maskFinal{frameNum}=mask;
    cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Detrmine ratio
    imFRETbg(~mask)=nan;
    imFRET=ndnanfilter(imFRETbg,fspecial('disk',3),'replicate');
    imCFPbg(~mask)=nan;
    imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');

    imRatio_raw{frameNum}=imFRET./imCFP;
    imRatio_raw{frameNum}(~mask)=nan;
   
    
    
   imRubybg(~mask)=nan;
  %  im_mRuby_=ndnanfilter(imRubybg,fspecial('disk',3),'replicate');
   im_mRuby_raw{frameNum}=imRubybg;
    
   imcytobg(~mask)=nan;
    im_cyto_raw{frameNum}=imcytobg;
   
   
    %%%%%% Generate and write files for raw ratio and outlined objects
   
    imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,mask);
    
    imwrite(imFRETOutline{frameNum},[datadir,filesep,position,'_Outline_prelim.tif'],'WriteMode','append','Compression','none');
end
%%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay
bleach_1=nanmean(vect(imRatio_raw{1}));
%bleach_raw = {}; bleach_raw_mRuby = {}; 
for frameNum=1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
  bleach_raw_mRuby(frameNum)= nanmean(vect(im_mRuby_raw{frameNum}));
  bleach_raw_cyto(frameNum)= nanmean(vect(im_cyto_raw{frameNum}));
end
save([datadir,filesep,position,'_RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw','imFRETOutline','-v7.3','im_mRuby_raw', 'im_cyto_raw'); 

save([datadir,filesep,position,'_Bleach_raw.mat'],'bleach_raw','bleach_raw_mRuby', 'bleach_raw_cyto');  %'CFP_raw','FRET_raw'); %
end  