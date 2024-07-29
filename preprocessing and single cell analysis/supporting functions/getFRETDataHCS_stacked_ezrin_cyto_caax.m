function getFRETDataHCS_stacked_ezrin_cyto_caax(cellNum,rawdir,datadir, threshold)

if ~exist([datadir,filesep,'RatioData_raw.mat'])
%%%%%% Set up

%load([rawdir,filesep,'alignment parameters pX pY.mat'],'pX','pY');

%%%%%% Call background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings

% cytobg_raw=double(imread([rawdir,filesep,'cyto_bg.tif']));
% caaxbg_raw=double(imread([rawdir,filesep,'caax_bg.tif']));
% mRubybg_raw=double(imread([rawdir,filesep,'mRuby_bg.tif']));

% bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
% bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
% CFPbg=bg2(:,:,1);
% FRETbg=bg2(:,:,2);
% 
%  [m,n]=size(mRubybg_raw);
%  bg3(:,:,1)=mRubybg_raw;bg3(:,:,2)=zeros(m,n);
%  bg4=dualviewAlignFromFittedSurface(bg3,pX,pY,binning);
%  mRubybg=bg4(:,:,1);

cellPath=strcat('cell_',num2str(cellNum));
cellFiles=getFilenames([rawdir],'.tif');

%  caax_stack=double(readTIFFstack([rawdir,filesep,cellFiles{1}]));
 cyto_stack=double(readTIFFstack([rawdir,filesep,cellFiles{1}]));
 mRuby_stack=double(readTIFFstack([rawdir,filesep,cellFiles{2}]));
 
  
%%%%%% Loop through frames
imRatio_raw={};maskFinal={};cellCoors={};   im_mRuby_raw={};
for frameNum=1:size(mRuby_stack,3)
    disp(num2str(frameNum));
    imcyto_raw=cyto_stack(:,:,frameNum);
%     imcaax_raw=caax_stack(:,:,frameNum);
    imRuby_raw=mRuby_stack(:,:,frameNum);

    %%%%%% Align CFP/FRET images
%     imstack(:,:,1)=imCFP_raw; imstack(:,:,2)=imFRET_raw;
%     imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
%     imCFP_raw=imaligned(:,:,1);
%     imFRET_raw=imaligned(:,:,2);
    
%       imstack2(:,:,1)=imRuby_raw; imstack2(:,:,2)=zeros(m,n);
%       imaligned2=dualviewAlignFromFittedSurface(imstack2,pX,pY,1);
%       imRuby_raw=imaligned2(:,:,1);
    %%%%%% Background-subtract CFP/FRET images
%     bgmask_caax=getBGMask(imcaax_raw);
%     bgmask_MRuby=getBGMask(imRuby_raw);
%      bgmask_cyto=getBGMask(imcyto_raw);% -- having issues here with
% 
%  imcaaxbg=subBG(imcaax_raw,bgmask_caax,caaxbg_raw);
%  imcytobg=subBG(imcyto_raw,bgmask_cyto,cytobg_raw);
%  imRubybg=subBG(imRuby_raw,bgmask_MRuby,mRubybg_raw);
    %%%%%% Get mask from raw FRET image
    
    s= cellNum; 
   [mask cellCoorsTemp]=getCellMaskCyto_2_stacked((2*imRuby_raw),1000, frameNum,s, threshold); % see here if its better taking away the imCFPbg
   
    maskFinal{frameNum}=mask;
   cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Detrmine ratio
%     imRubybg(~mask)=nan;
%      im_mRubyTEMP=ndnanfilter(imRubybg,fspecial('disk',1),'replicate');
%       im_mRubyTEMP(~mask)=nan;
%      im_mRuby_raw{frameNum}=im_mRubyTEMP;
imRuby_raw(~mask)=nan;
     im_mRubyTEMP=ndnanfilter(imRuby_raw,fspecial('disk',1),'replicate');
      im_mRubyTEMP(~mask)=nan;
     im_mRuby_raw{frameNum}=im_mRubyTEMP;
% 
%       imcytobg(~mask)=nan;
%    im_cytoTEMP =ndnanfilter(imcytobg,fspecial('disk',1),'replicate');
%       im_cytoTEMP(~mask)=nan; 
%     im_cyto_raw{frameNum}=im_cytoTEMP;
    
     imcyto_raw(~mask)=nan;
   im_cytoTEMP =ndnanfilter(imcyto_raw,fspecial('disk',1),'replicate');
      im_cytoTEMP(~mask)=nan; 
    im_cyto_raw{frameNum}=im_cytoTEMP;
    
%     imcaaxbg(~mask)=nan;
%     im_caax_raw{frameNum}=imcaaxbg;
   
   
%     imCAAXOutline{frameNum}=DrawMaskOutline(imcaax_raw,mask);
 immRuby_outline{frameNum}=DrawMaskOutline(imRuby_raw,mask);
     

  % imwrite(imCAAXOutline{frameNum},[datadir,filesep,'Outline_prelim.tif'],'WriteMode','append','Compression','none');
      imwrite(immRuby_outline{frameNum},[datadir,filesep,'Outline_mRuby_prelim.tif'],'WriteMode','append','Compression','none');
    
end
%%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay

for frameNum= 1:length(im_mRuby_raw)
%    bleach_raw_caax(frameNum)=nanmean(vect(im_caax_raw{frameNum}));
   bleach_raw_mRuby(frameNum)= nanmean(vect(im_mRuby_raw{frameNum}));
   bleach_raw_cyto(frameNum)= nanmean(vect(im_cyto_raw{frameNum}));
end
save([datadir,filesep,'RatioData_raw.mat'],'maskFinal','cellCoors','im_mRuby_raw', 'im_cyto_raw','-v7.3'); %'im_caax_raw', 'imCAAXOutline'
save([datadir,filesep,'Bleach_raw.mat'],'bleach_raw_mRuby', 'bleach_raw_cyto'); %'bleach_raw_caax',
end
