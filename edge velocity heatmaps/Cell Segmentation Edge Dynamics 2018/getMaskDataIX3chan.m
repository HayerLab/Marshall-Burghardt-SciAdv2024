function getMaskDataIX(position,bgdir,rawdir,datadir,jitter)

% image processing for FRET data analysis 
% row=3;col=3;site=2;
% Input: raw images, background images, and alignment parameters for
% images acquired in different channels. 
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
%   ndnanfilter
%
% Arnold, 180622

if ~exist([datadir,filesep,position,'_RatioData_raw.mat'])

load([bgdir,filesep,'alignment parameters pX pY.mat'],'pX','pY');

%%%%%% Call background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
TRITCbg_raw=double(imread([bgdir,filesep,'AVG_bgTRITC.tif']));
CFPbg_raw=double(imread([bgdir,filesep,'AVG_bgCFP.tif']));
YFPbg_raw=double(imread([bgdir,filesep,'AVG_bgYFP.tif']));
%align CFP to TRITC, crop YFP background
bg1(:,:,1)=TRITCbg_raw; bg1(:,:,2)=CFPbg_raw;
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
TRITCbg=bg2(:,:,1);
CFPbg=bg2(:,:,2);

bg3(:,:,1)=YFPbg_raw; bg3(:,:,2)=TRITCbg_raw;
bg4=dualviewAlignFromFittedSurface(bg3,pX,pY,binning);
YFPbg=bg4(:,:,1);

TRITC_files=dir([rawdir,filesep,position,'_TRITC_*']);

%%%%%% Loop through frames
imRatio_raw={};maskFinal={};cellCoors={};
for frameNum=1:length(TRITC_files)
    disp([position,'__',num2str(frameNum)]);
    imTRITC_raw=double(imread([rawdir,filesep,position,'_TRITC_',num2str(frameNum),'.tif']));
    imCFP_raw=double(imread([rawdir,filesep,position,'_CFP_',num2str(frameNum),'.tif']));
    imYFP_raw=double(imread([rawdir,filesep,position,'_YFP_',num2str(frameNum),'.tif']));
    
    %%%%%% Align TRITC/CFP images
    imstack(:,:,1)=imTRITC_raw; imstack(:,:,2)=imCFP_raw;
    imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
    imTRITC_al=imaligned(:,:,1);
    imCFP_al=imaligned(:,:,2);
    
    % Generate cropped YFP frame as well
    imstack2(:,:,1)=imYFP_raw;
    imstack2(:,:,2)=imTRITC_raw;
    imaligned2=dualviewAlignFromFittedSurface(imstack2,pX,pY,binning);
    imYFP_al=uint16(imaligned2(:,:,1));
        
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imTRITC_al+imCFP_al);
    imTRITCbg=subBG(imTRITC_al,bgmask,TRITCbg);
    imCFPbg=subBG(imCFP_al,bgmask,CFPbg);
    imYFPbg=subBG(imYFP_al,bgmask,YFPbg);
    %%%%%% Get mask from raw TRITC image
    [mask cellCoorsTemp]=getCellMask(imTRITC_al,2000);
    maskFinal{frameNum}=mask;
    cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Determine ratio
    imCFPbg(~mask)=nan;
    imCFP=ndnanfilter(imCFPbg,fspecial('disk',3),'replicate');
    imTRITCbg(~mask)=nan;
    imTRITC=ndnanfilter(imTRITCbg,fspecial('disk',3),'replicate');
    imRatioTemp=imCFP_al./imTRITC_al;
    imRatioTemp(~mask)=nan;
    imRatio_raw{frameNum}=imRatioTemp;
    
    imYFPbg(~mask)=nan;
    imYFP=ndnanfilter(imYFPbg,fspecial('disk',2),'replicate');
    imYFPcrop{frameNum}=imYFP;
    
    %%%%%% Determine scaling for representation
    if frameNum==1
       colorRange = [round(prctile(imRatioTemp(:),3),1),round(prctile(imRatioTemp(:),97),1)];
    end
    %%%%%% Generate and write files for raw ratio and outlined objects
    tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imTRITCOutline{frameNum}=DrawMaskOutline(imTRITC_al,mask);
    imwrite(imTRITCOutline{frameNum},[datadir,filesep,position,'_TRITC.tif'],'WriteMode','append','Compression','none');
    imwrite(tempRATIO,[datadir,filesep,position,'_RATIO_raw_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
%    
end
%%%%%% Determine Bleaching over the course of the experiment (fluorescence
%%%%%% decay in masked areas.
bleach_1=nanmean(vect(imRatio_raw{1}));
for frameNum=1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
end
save([datadir,filesep,position,'_MaskingData.mat'],'maskFinal','cellCoors','imRatio_raw','imYFPcrop','-v7.3');
save([datadir,filesep,position,'_Bleach_raw.mat'],'bleach_raw');
end



