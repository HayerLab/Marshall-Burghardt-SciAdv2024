function getFRETDataHCS_stacked_3chan(cellNum,rawdir,datadir, threshold, pX, pY)

if ~exist([datadir,filesep,'RatioData_raw.mat'])
%%%%%% Set up



%%%%%% Call background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
CFPbg_raw=double(imread([rawdir,filesep,'CFP_bg.tif']));
FRETbg_raw=double(imread([rawdir,filesep,'FRET_bg.tif']));
mRubybg_raw=double(imread([rawdir,filesep,'mRuby_bg.tif']));

bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);
% 
 [m,n]=size(mRubybg_raw);
 bg3(:,:,1)=mRubybg_raw;bg3(:,:,2)=zeros(m,n);
 bg4=dualviewAlignFromFittedSurface(bg3,pX,pY,binning);
 mRubybg=bg4(:,:,1);

cellPath=strcat('cell_',num2str(cellNum));
cellFiles=getFilenames([rawdir],'.tif');

 CFP_stack=double(readTIFFstack([rawdir,filesep,cellFiles{1}]));
 FRET_stack=double(readTIFFstack([rawdir,filesep,cellFiles{3}]));
 mRuby_stack=double(readTIFFstack([rawdir,filesep,cellFiles{5}]));
 
  
%%%%%% Loop through frames
imRatio_raw={};maskFinal={};cellCoors={};   im_mRuby_raw={};
for frameNum=1:size(CFP_stack,3)
    disp(num2str(frameNum));
    imCFP_raw=CFP_stack(:,:,frameNum);
    imFRET_raw=FRET_stack(:,:,frameNum);
    imRuby_raw=mRuby_stack(:,:,frameNum);

    %%%%%% Align CFP/FRET images
    imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
    imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
    imCFP_raw=imaligned(:,:,1);
    imFRET_raw=imaligned(:,:,2);
    
      imstack2(:,:,1)=imRuby_raw; imstack2(:,:,2)=zeros(m,n);
      imaligned2=dualviewAlignFromFittedSurface(imstack2,pX,pY,1);
      imRuby_raw=imaligned2(:,:,1);
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imCFP_raw+imFRET_raw);
    bgmask_MRuby=getBGMask(imRuby_raw); 

 imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
 imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
 imRubybg=subBG(imRuby_raw,bgmask_MRuby,mRubybg);
    %%%%%% Get mask from raw FRET image
    
    s= cellNum; 
   [mask cellCoorsTemp]=getCellMaskCyto_2_stacked((2*imFRETbg),1000, frameNum,s, threshold); 
   
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
    
     imRubybg(~mask)=nan;
   
     %if you desire filtering of mRuby channel 
     %im_mRuby_raw=ndnanfilter(imRubybg,fspecial('disk',3),'replicate');
     im_mRuby_raw{frameNum}=imRubybg;


   
    imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,mask);
    immRuby_outline{frameNum}=DrawMaskOutline(imRuby_raw,mask); 

   imwrite(imFRETOutline{frameNum},[datadir,filesep,'Outline_prelim.tif'],'WriteMode','append','Compression','none');
   %imwrite(immRuby_outline{frameNum},[datadir,filesep,'Outline_mRuby_prelim.tif'],'WriteMode','append','Compression','none');
    
end
%%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay

for frameNum= 1:length(imRatio_raw)
   bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
   bleach_raw_mRuby(frameNum)= nanmean(vect(im_mRuby_raw{frameNum}));
end
save([datadir,filesep,'RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw','im_mRuby_raw','immRuby_outline', 'imFRETOutline','-v7.3'); 
save([datadir,filesep,'Bleach_raw.mat'],'bleach_raw','bleach_raw_mRuby');
end
