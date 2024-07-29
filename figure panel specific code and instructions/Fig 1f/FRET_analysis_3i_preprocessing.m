% Preprocessing for ratiometric FRET analysis. 
% 
% - Generation of averaged images of frame 3 of all CFP and FRET channels 
%   time-lapse sequences, to be used for image registration.
%   
% - Determination of parameters for image registration between CFP and FRET
%   channels, developed by Sean Collins and described in Yang et al. 2016, 
%   NCB (PMID:26689677) 
%
% Requires the following not built-in functions:
% 
%   getFilenames.m
%   boolREgExp.m
%   readTIFFstack.m
%   dualviewComputeAlignmentFromGridImages.m
%   dualview2stack.m
%   getOffsetBetweenImages.m
%   dualviewGetXMatGridFor2ndOrderSurface.m
%   maxMatrixCoord.m

%% Preparation
clear 
clc 
close all
warning off
root = 'Y:\Rho paper data\MainFigures\Fig1\Fig 1f\170816_3i\Thrombin';
rawdir=[root,filesep,'raw']; % folder where raw image TIFF stacks are stored
bgdir=[root,filesep,'background'];

%% Generate averaged images of third frame for alignment
CFP_stack=[];FRET_stack=[];

files=getFilenames(rawdir);
files=files(boolRegExp(files,'.tiff'));
CFPfiles=files(boolRegExp(files,'_CFP'));
FRETfiles=files(boolRegExp(files,'_FRET'));

CFP_stack=[];FRET_stack=[];
for numfiles=1:length(CFPfiles)
   disp(CFPfiles(numfiles));
   CFP_temp=readTIFFstack([rawdir,filesep,CFPfiles{numfiles}]);
   FRET_temp=readTIFFstack([rawdir,filesep,FRETfiles{numfiles}]);
   CFP_stack(:,:,numfiles)=CFP_temp(:,:,3);
   FRET_stack(:,:,numfiles)=FRET_temp(:,:,3);
end

CFP_AV=uint16(mean(CFP_stack,3));
FRET_AV=uint16(mean(FRET_stack,3));
imwrite(CFP_AV,[bgdir,filesep,'AVG_rawdata_CFP.tif'],'TIFF','Compression','None');
imwrite(FRET_AV,[bgdir,filesep,'AVG_rawdata_FRET.tif'],'TIFF','Compression','None');

%% Determine parameters for image registration
imCFP_raw=imread([bgdir,filesep,'AVG_rawdata_CFP.tif']);
imFRET_raw=imread([bgdir,filesep,'AVG_rawdata_FRET.tif']);
alignStack(:,:,1)=imCFP_raw; alignStack(:,:,2)=imFRET_raw;
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

save([bgdir,filesep,'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');
disp('done!');
