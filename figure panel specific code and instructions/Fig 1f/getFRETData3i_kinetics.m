function getFRETData3i_kinetics(position,bgdir,rawdir,datadir,jitter)

% image processing for FRET data analysis 

% Input: CFP/FRET images, background images, and alignment parameters for
% CFP/FRET images. 
%
% Generates first a mask for background subtraction, subtracts background,
% and calculates a refined mask before taking a ra
% tio between background-
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
%   imRatio_raw Cell array with masked raw ratio values, before bleaching
%               correction
%
%   imFRETOutline.tif 
%               RGB tif stack of FRET channel images with outlined masked
%               areas
%               %
%
% Used non built-in functions: 
%
%   dualviewAlignFromFittedSurface.m
%   getBGMask.m
%   subBG.m
%   getCellMask_kinetics.m
%   DrawMaskOutline.m
%   ratio2RGB.m
%   vect.m
%   readTIFFstack.m
%   ndnanfilter.m

% Arnold Hayer

load([bgdir,filesep,'alignment parameters pX pY.mat'],'pX','pY');

binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
CFPbg_raw=double(imread([bgdir,filesep,'AVG_bg_CFP.tif']));
FRETbg_raw=double(imread([bgdir,filesep,'AVG_bg_FRET.tif']));
bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
CFPbg=bg2(:,:,1);
FRETbg=bg2(:,:,2);

files=getFilenames(rawdir);
files=files(boolRegExp(files,'.tiff'));
CFPfiles=files(boolRegExp(files,'_CFP'));
FRETfiles=files(boolRegExp(files,'_FRET'));
%REDfiles=files(boolRegExp(files,'_TRITC'));
%REDfiles=files(boolRegExp(files,'_Cherry'));

fileID=CFPfiles{position};
fileID=fileID(1:(end-5)); % removes '.tiff' from string 

if ~exist([datadir,filesep,fileID,'_RatioData_raw.mat'])
    disp(fileID);
    
    CFP_stack=double(readTIFFstack([rawdir,filesep,CFPfiles{position}]));
    FRET_stack=double(readTIFFstack([rawdir,filesep,FRETfiles{position}]));
    %RED_stack=double(readTIFFstack([rawdir,filesep,REDfiles{position}]));
        
    %%%%%% Loop through frames
    imRatio_raw={};maskFinal={};cellCoors={};
    for frameNum=1:size(CFP_stack,3)
        disp([num2str(position),'__',num2str(frameNum)]);
           
        %%%%%% Align CFP/FRET images
        imstack(:,:,1)=CFP_stack(:,:,frameNum);
        imstack(:,:,2)=FRET_stack(:,:,frameNum);
        imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,binning);
        imCFP_raw=imaligned(:,:,1);
        imFRET_raw=imaligned(:,:,2);
        
        %% Generate cropped RED frame as well
        %imstack2(:,:,1)=RED_stack(:,:,frameNum);
        %imstack2(:,:,2)=FRET_stack(:,:,frameNum);
        %imaligned2=dualviewAlignFromFittedSurface(imstack2,pX,pY,binning);
        %imRED_raw=uint16(imaligned2(:,:,1));
        
        %%%%%% Background-subtract CFP/FRET images
        bgmask=getBGMask(imCFP_raw+imFRET_raw);
        imCFPbg=subBG(imCFP_raw,bgmask,CFPbg);
        imFRETbg=subBG(imFRET_raw,bgmask,FRETbg);
        %%%%%% Get mask from raw FRET image
        [mask cellCoorsTemp]=getCellMask_kinetics(imFRET_raw+imCFP_raw,4000);
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
        %%%%%% Determine scaling for representation
        if frameNum==1
           colorRange = [round(prctile(imRatioTemp(:),3),1),round(prctile(imRatioTemp(:),97),1)];
        end
        %%%%%% Generate datastructure with outlined objects based on mask
        imFRETOutline{frameNum}=DrawMaskOutline(imFRET_raw,mask);
       
    end
    %%%%%% Compute overall decay kinetics   
    for frameNum=1:length(imRatio_raw)
       bleach_raw(frameNum)=nanmean(vect(imRatio_raw{frameNum}));
    end

    save([datadir,filesep,fileID,'_RatioData_raw.mat'],'maskFinal','cellCoors','imRatio_raw','imFRETOutline','-v7.3');
    save([datadir,filesep,fileID,'_Bleach_raw.mat'],'bleach_raw');
end


