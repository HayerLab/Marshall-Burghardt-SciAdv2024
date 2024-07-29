clear;clc;close all; 
root = 'H:\180620_IXM\180620-I';
rawdir=[root,filesep,'raw'];
bgdir=[root,filesep,'background'];
load([bgdir,filesep,'alignment parameters pX pY.mat']);
datadir=[rawdir,filesep,'data_180622'];
position='5_8_1';

%% Load and process background images
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
TRITCbg_raw=double(imread([bgdir,filesep,'AVG_bgTRITC.tif']));
GFPbg_raw=double(imread([bgdir,filesep,'AVG_bgGFP.tif']));
bg1(:,:,1)=TRITCbg_raw; bg1(:,:,2)=GFPbg_raw;
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
TRITCbg=bg2(:,:,1);
GFPbg=bg2(:,:,2);
TRITC_files=dir([rawdir,filesep,position,'_TRITC_*']);

%%
imRatio_raw={};maskFinal={};cellCoors={};
    frameNum=1%:length(TRITC_files)
    disp([position,'__',num2str(frameNum)]);
    imTRITC_raw=double(imread([rawdir,filesep,position,'_TRITC_',num2str(frameNum),'.tif']));
    imGFP_raw=double(imread([rawdir,filesep,position,'_GFP_',num2str(frameNum),'.tif']));

    %%%%%% Align TRITC/GFP images
    imstack(:,:,1)=imTRITC_raw; imstack(:,:,2)=imGFP_raw;
    imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
    imTRITC_al=imaligned(:,:,1);
    imGFP_al=imaligned(:,:,2);
    %%%%%% Background-subtract CFP/FRET images
    bgmask=getBGMask(imTRITC_al+imGFP_al);
    imTRITCbg=subBG(imTRITC_al,bgmask,TRITCbg);
    imGFPbg=subBG(imGFP_al,bgmask,GFPbg);
    %%%%%% Get mask from raw FRET image
    [mask cellCoorsTemp]=getCellMask(imTRITC_al,2000);
    maskFinal{frameNum}=mask;
    cellCoors{frameNum}=cellCoorsTemp;
    %%%%%% Determine ratio
    imGFPbg(~mask)=nan;
    imGFP=ndnanfilter(imGFPbg,fspecial('disk',3),'replicate');
    imTRITCbg(~mask)=nan;
    imTRITC=ndnanfilter(imTRITCbg,fspecial('disk',3),'replicate');
    imRatioTemp=imGFP_al./imTRITC_al;
    imRatioTemp(~mask)=nan;
    imRatio_raw{frameNum}=imRatioTemp;
    %%%%%% Determine scaling for representation
    if frameNum==1
       colorRange = [round(prctile(imRatioTemp(:),3),1),round(prctile(imRatioTemp(:),97),1)];
    end
    %%%%%% Generate and write files for raw ratio and outlined objects
    tempRATIO=ratio2RGB(imRatioTemp,colorRange);
    imTRITCOutline{frameNum}=DrawMaskOutline(imTRITC_al,mask);
%     imRFPOutline=DrawMaskOutlineChy(imRFP_raw,mask);
    %imwrite(imFRETOutline,[datadir,filesep,position,'_FRET.tif'],'WriteMode','append','Compression','none');
    %imwrite(tempRATIO,[datadir,filesep,position,'_RATIO_raw_',num2str(colorRange(1)),'_',num2str(colorRange(2)),'.tif'],'WriteMode','append','Compression','none');
%     imwrite(imRFPOutline,[datadir,filesep,shot,'_RFP.tif'],'WriteMode','append','Compression','none');

%% Display alignment
figure;
subplot(1,2,1);showImagesMergeChannels(imstack(:,:,1),imstack(:,:,2));
subplot(1,2,2);showImagesMergeChannels(imTRITC_al,imGFP_al);

%% Display segmentation
figure;
testimage=[imTRITCOutline{frameNum} tempRATIO];
imshow(testimage);

%% background mask
figure; imagesc(bgmask);


