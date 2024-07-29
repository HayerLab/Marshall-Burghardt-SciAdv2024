%% Purpose: To determine alignment parameteters between two fluorescence channels
%      that are not perfectly aligned due to chromatic aberration or imperfections
%      of alingment of optical components in the light path. 
% 
%  Input: 
%      Timelapse sequences of images acquired using an IX micro
%      fluorescence microscope with filenames in the format
%      row_col_site_channel1_timepoint.tif and row_col_site_channel2_timepoint.tif
% 
%  Output: 
%      image alignment parameters, stored in 'alignment parameters pX pY.mat'
% 
%  Required functions: 
%      dualviewAlignFromFittedSurface.m
%      dualviewComputeAlignmentFromGridImages.m
%      dualviewGetXMatGridFor2ndOrderSurface.m
%      getOffsetBetweenImages.m
%      dualviewAlignFromFittedSurface.m
%      showImagesWithLinkedAxes.m
%      dualview2stack.m
%      maxMatrixCoord.m     
% Arnold Hayer, 22 June 2018

%% Working folder
clear 
clc 
close all
%warning off
root = 'H:\180620_IXM\180620-I';
rawdir=[root,filesep,'raw'];
bgdir=[root,filesep,'background'];
chan1='TRITC';
chan2='GFP';

%% Generate averaged image of third (arbitrary) frame for alignment
chan1_stack=[];chan2_stack=[];i=0;
for row=5%2:3
    for col=8:10
        for site=1:4
            i=i+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            chan1_temp=imread([rawdir,filesep,shot,'_',chan1,'_3.tif']);
            chan2_temp=imread([rawdir,filesep,shot,'_',chan2,'_3.tif']);
            chan1_stack(:,:,i)=chan1_temp;
            chan2_stack(:,:,i)=chan2_temp;
        end 
    end
end
chan1_AV=uint16(mean(chan1_stack,3));
chan2_AV=uint16(mean(chan2_stack,3));
imwrite(chan1_AV,[bgdir,filesep,'AVG_rawdata_',chan1,'.tif'],'TIFF','Compression','None');
imwrite(chan2_AV,[bgdir,filesep,'AVG_rawdata_',chan2,'.tif'],'TIFF','Compression','None');
% 
%% Determine image alignment using Sean's method
im1_raw=imread([bgdir,filesep,'AVG_rawdata_',chan1,'.tif']);
im2_raw=imread([bgdir,filesep,'AVG_rawdata_',chan2,'.tif']);
alignStack(:,:,1)=im1_raw; alignStack(:,:,2)=im2_raw;
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

save([bgdir,filesep,'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');

%% Test parameters
load([bgdir,filesep,'alignment parameters pX pY.mat'],'pX','pY');
frameNum=1;
%%%%%%%%%%%%%% Specify row, col, site to test:  
shot='7_11_3';
%%%%%%%%%%%%%%
chan1_files=dir([rawdir,filesep,shot,'_',chan1,'_*']);
chan2_files=dir([rawdir,filesep,shot,'_',chan2,'_*']);
im1_raw=double(imread([rawdir,filesep,shot,'_',chan1,'_',num2str(frameNum),'.tif']));
im2_raw=double(imread([rawdir,filesep,shot,'_',chan2,'_',num2str(frameNum),'.tif']));
imstack(:,:,1)=im1_raw; imstack(:,:,2)=im2_raw;
imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);
subplot(1,2,1), showImagesMergeChannels(im1_raw,im2_raw);
subplot(1,2,2), showImagesMergeChannels(imaligned(:,:,1),imaligned(:,:,2));

%% Load and process background image (an average of multiple background images)
binning=1; % relevant if alingment images and data images were acquired using distinct binning settings
im1bg_raw=double(imread([bgdir,filesep,'AVG_bg',chan1,'.tif']));
im2bg_raw=double(imread([bgdir,filesep,'AVG_bg',chan2,'.tif']));
bg1(:,:,1)=im1bg_raw; bg1(:,:,2)=im2bg_raw;
bg2=dualviewAlignFromFittedSurface(bg1,pX,pY,binning);
im1bg=bg2(:,:,1);
im2bg=bg2(:,:,2);

%% Test parameters and alignment using a single frame
close all
%%%%%%%%%%%%%% Specify row,col,site here
shot= '7_11_1';
%%%%%%%%%%%%%%
frameNum=1;
im1_raw=double(imread([rawdir,filesep,shot,'_',chan1,'_',num2str(frameNum),'.tif']));
im2_raw=double(imread([rawdir,filesep,shot,'_',chan2,'_',num2str(frameNum),'.tif']));
imstack(:,:,1)=im1_raw; imstack(:,:,2)=im2_raw;
imaligned=dualviewAlignFromFittedSurface(imstack,pX,pY,1);

im1_raw=imaligned(:,:,1);
im2_raw=imaligned(:,:,2);

bgmask=getBGMask(im1_raw+im2_raw);
im1bg=subBG(im1_raw,bgmask,im1bg);
im2bg=subBG(im2_raw,bgmask,im2bg);

subplot(1,2,2);imagesc(im2bg+im1bg);
subplot(1,2,1);imagesc(bgmask);


