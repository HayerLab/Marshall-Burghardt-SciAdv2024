function getMaskData(position,rawdir,datadir)
%
% 
% Input: Time-series of fluorescenct images, background images (optional)
%
% Generates first (optional) a mask for background subtraction, subtracts background,
% and calculates a refined mask 
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
%   getBGMask
%   subBG
%   getCellMask
%   DrawMaskOutline
%   ratio2RGB
%   vect
%   readTIFFstack

% Arnold, 24 May 2015




%%%%%% Set up
% root = 'K:\data\150424_IX\RhoA2G-HT73-30s';
% rawdir=[root,filesep,'rawdata'];
% % bgdir='K:\data\150424_IX\background-bin1';
% %load([bgdir,filesep,'alignment parameters pX pY.mat'],'pX','pY');
% datadir=[root,filesep,'data_150527'];
% shot=[num2str(row),'_',num2str(col),'_',num2str(site)];


%%%%%% Call background images
% binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
% CFPbg_raw=double(imread([bgdir,filesep,'AVG_bg_CFP.tif']));
% FRETbg_raw=double(imread([bgdir,filesep,'AVG_bg_FRET.tif']));
% bg1(:,:,1)=CFPbg_raw; bg1(:,:,2)=FRETbg_raw;
% bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
% CFPbg=bg2(:,:,1);
% FRETbg=bg2(:,:,2);

% With or withoug background substraction
 %backgroundsubtraction == FALSE; -- not needed here

files=getFilenames(rawdir);
%files=files(boolRegExp(files,'.tif'));
%files=files(boolRegExp(files,'Wide_5'));
% YFPfiles=files(boolRegExp(files,'_YFP'));
%RFPfiles=files(boolRegExp(files,'_Cherry'));
%REDfiles=files(boolRegExp(files,'_TRITC'));
%REDfiles=files(boolRegExp(files,'_mCherry'));

fileID=files{position};
fileID=fileID(1:(end-4)); % removes '.tif' from string 

if ~exist([datadir,filesep,fileID,'_MaskData5.mat'])
 
    
    YFP_stack=double(readTIFFstack([rawdir]));  %,filesep,files{position}]));
%    RFP_stack=double(readTIFFstack([rawdir,filesep,RFPfiles{position}]));
%    RED_stack=double(readTIFFstack([rawdir,filesep,REDfiles{position}]));
        
    %%%%%% Loop through frames
    maskFinal={};cellCoors={};
    for frameNum=1:size(YFP_stack,3)
        disp([num2str(position),'__',num2str(frameNum)]);
     
      imYFP_raw=YFP_stack(:,:,frameNum);
       % imRFP_raw=RFP_stack(:,:,frameNum);
        
%         if backgroundsubstraction
%             %%%%%% Background-subtract 
%             bgmask=getBGMask(imYFP_raw);
%             imYFPbg=subBG(imYFP_raw,bgmask);
%         end
%         
        %%%%%% Get mask from raw image
        
       % bg_raw = double(imread(['F:\Seph\data\data_200610 - Trial 1 Rac Myosin\background\AVG_bgCFP.tif'])); 
      bg_raw = double(imread(['F:\Seph\data\data_200610 - Trial 1 Rac Myosin\background\AVG_bgCFP_cropped_20.tif']));
         bgmask_mRuby=getBGMask(imYFP_raw);
    
    imRubybg=subBG(imYFP_raw,bgmask_mRuby,bg_raw);
        
       [mask cellCoorsTemp]=getCellMaskCyto_2(imRubybg,2000); %turned this into get cell mask cyto
        maskFinal{frameNum}=mask;
        cellCoors{frameNum}=cellCoorsTemp;
       
        
        imYFPOutline{frameNum}=DrawMaskOutline(imYFP_raw,mask);
        %imRFPOutline{frameNum}=DrawMaskOutline(imRFP_raw,mask);
        %stitched=[imYFPOutline{frameNum} imRFPOutline{frameNum}];
        imwrite(imYFPOutline{frameNum},[datadir,filesep,fileID,'_masked-using 2 function.tif'],'WriteMode','append','Compression','none');
        %imRFPOutline=DrawMaskOutlineChy(imRFP_raw,mask);
        %imwrite(imFRETOutline,[datadir,filesep,position,'_FRET.tif'],'WriteMode','append','Compression','none');
      
    end
    %%%%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay
  
%     fileID=CFPfiles{position};
%     fileID=fileID(1:(end-5)); % removes '.tiff' from string 
    save([datadir,filesep,fileID,'_MaskData2.mat'],'maskFinal','cellCoors','imYFPOutline','-v7.3');
    

end
end


function [maskFinal cellCoors] = getCellMask(image,  minCellSize)
% getMask returns a binary mask based on image and minCellSize
% Uses a histogram-based thresholding
% Stores centroid and area of detected objects
% Arnold Hayer 150520

% if nargin<3
%     maxCellSize=16000;
% end
if nargin<2
    minCellSize=4000;
end
% Subtract scaled blurred image to sharpen outline
imSmooth1=imfilter(image,fspecial('disk',1),'replicate');
image_lpf=imfilter(imSmooth1,fspecial('disk',20),'symmetric');
imSmooth2=imSmooth1-image_lpf*0.6; % arbitrary scaling factor

% Histogram-based threshold determination
% find range of pixel values
PixMinMax=double([0 round(prctile(imSmooth2(:),99))]);
IntStep=ceil((PixMinMax(2)-PixMinMax(1))/300);  
[f,xi]=ksdensity(imSmooth2(:),PixMinMax(1):IntStep:PixMinMax(2));

%figure; plot(xi,f); hold on;
    
[pks,locs]=findpeaks(f);
%filter peaks greater than 0.005
log_pks=pks>0.0005; 

pks=pks(log_pks); locs=locs(log_pks);

%plot(xi(locs),pks,'k^','markerfacecolor',[1 0 0]);hold off;
x_bgMax=xi(locs(1)); % picks the first peak (x-value of first peak)

[~,ind]=find(f>(0.01*pks(1)),1); % returns the first value of f greater than 0.5% of its max
x_1pct=xi(ind); % returns the corresponding intensity value
bgWidth=x_bgMax-x_1pct; % estimates the width of the background peak
threshSeg=(x_bgMax+bgWidth); % multiplied by 3 because of good fg/bg separationthreshold level based on the background peak
%figure, showImagesWithLinkedAxes({imSum,imerode(imSumSmooth>threshSeg,strel('disk',1))});

%mask_init=imerode(imSmooth2>threshSeg,strel('disk',2));
mask_init=imSmooth2>threshSeg;
%mask_holes_filled=imfill(mask_init,'holes');
% With or without filling holes
%mask=bwareaopen(mask_holes_filled,minCellSize);


% Display result
%subplot(2,2,1); imshow(mat2gray(image));
%subplot(2,2,2);imshow(mat2gray(mask));

%maskinv=imcomplement(mask);
seD = strel('diamond',1);% was 2
mask = imerode(mask_init,seD);
mask_holes_filled=imfill(mask,'holes');

% maskBlobs=bwareaopen(mask_holes_filled,maxCellSize);
% mask_holes_filled(maskBlobs)=0;
maskFinal=bwareaopen(mask_holes_filled,minCellSize);

% BWoutline = bwperim(maskFinal);
% Segout = image;
% Segout(BWoutline) = max(image(:));
% figure, imagesc(Segout);

% get centroid coordinates
objects=regionprops(maskFinal,'BoundingBox','Area','Centroid','PixelIdxList','PixelList');
objects=objects(arrayfun(@(x) x.Area>minCellSize,objects));
%objects=objects(arrayfun(@(x) x.Area<maxCellSize,objects));

cellCoors(:,1)=vect(arrayfun(@(x) x.Centroid(1),objects));
cellCoors(:,2)=vect(arrayfun(@(x) x.Centroid(2),objects));
cellCoors(:,3)=vect(arrayfun(@(x) x.Area,objects));
end
