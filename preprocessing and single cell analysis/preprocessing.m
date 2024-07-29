%% master script
% intialization
%ND2 to TIFF
% tiff stacks 
% background alignment 
% create RatioData_raw and bleaching correction curves 

% to run properly, store ND2 files in folder named 'ND2 files' with
% run1.nd2, run2.nd2 etc. This puts them in alphabetical
% order. Background ND2 file must be stored in fodler named 'bg ND2 files'

% before running any code here, refer to the text file 'File Structure
% instructions'. All code here is built for the specific naming system
% outlined there. 

% Seph Marshall Feb 2023

%% setup and initialization. 
clc; 
clear; 
addpath('preprocessing and single cell analysis'); 

%manually add root path where your raw data is stored 
root= 'F:\example_dataset_240522'; 

% if you want to create Tiff Stacks. Default set to yes, 1. 
makeTiffStacks = 1; 

% initialize number of channels along with name 
channels = {'FRET';'CFP'; 'mRuby';};
channel1 = channels'; 

channel = [ "FRET"; "CFP"; "mRuby"];


%input number of frames in the time series. 1 for single frame  
frames = 60; 

bgdir = [root, filesep, 'background']; 
if  ~exist(bgdir)
    mkdir(bgdir)
end 

rawdir = [root, filesep, 'raw']; 
if  ~exist(rawdir)
    mkdir(rawdir)
end 

datadir = [root, filesep, 'data']; 
if  ~exist(datadir)
    mkdir(datadir)
end 

bgND2dir = [root, filesep, 'ND2 files background']; 
if  ~exist(bgND2dir)
    mkdir(bgND2dir)
end 

ND2dir = [root, filesep, 'ND files']; 
if  ~exist(bgND2dir)
    mkdir(bgND2dir)
end 


if makeTiffStacks ==1 
    TIFFdir = [root, filesep, 'tiff_stacks']; 
if  ~exist(TIFFdir)
    mkdir(TIFFdir)
end 
end 
%% background ND2 to TIFF
% converts background ND2 files to tiff format. Here, ensure that for the
% correct number of channels, you have named them as desired and in the correct 
% order that they were imaged. 

bgND2files = getFilenames([bgND2dir],'.nd2');
bgfilepath = [bgND2dir, filesep, bgND2files{1}]; 
counter =0; 


finfo = nd2finfo(bgfilepath);
num_sites = finfo.img_seq_count; 

for row = 1
    for col = 1
        for site =1:num_sites
            
            
          for timept = 1:1
               
            shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
            
        if size (channels, 1) == 1
             [mCit]= nd2read_hayer_1chan(bgfilepath,finfo, counter,1); 
             imwrite(mCit,[bgdir,filesep,shot,'_mCit_',num2str(timept),'.tif'],'TIFF','Compression','None');
        end 
            
       if size(channels, 1) ==2
         [FRET, CFP]= nd2read_hayer_2chan(bgfilepath,finfo, counter,1); 
         imwrite(FRET,[bgdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
         imwrite(CFP,[bgdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
       end 
       
       if size (channels,1) == 3
              [FRET, CFP, mRuby]= nd2read_hayer(bgfilepath,finfo, counter,1); 
              imwrite(FRET,[bgdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(CFP,[bgdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[bgdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
  
       end 
            
       if size (channels,1) == 4
              [FRET,CFP,mRuby,cyto]= nd2read_hayer(bgfilepath,finfo, counter,1); 
        imwrite(CFP,[bgdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(FRET,[bgdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[bgdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
         imwrite(cyto,[bgdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
            
      if size (channels,1) == 6
             [mRuby, blank, cyto, blank2, FRET, CFP ] = nd2read_hayer_6chan(filepath,finfo,counter,1);
    imwrite(CFP,[bgdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(FRET,[bgdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[bgdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(cyto,[bgdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
      end 
      %
      
      
    
      
      counter = counter+1;
      
            end   
         end 
        end 
    end 
end 

%% ND2 file conversion to tiff for raw data 
% strongly recommend maunally checking that the raw files have been saved
% with the correct name and order. Occasionaly, files are saved in a
% different order from chich they were imaged. The only thing that needs to
% be changed is the order of the names in the output of the nd2read_hayer_
% functions. 

ND2files = getFilenames([ND2dir],'.nd2');


 for i = 1:(size(ND2files,1))
filepath = [ND2dir, filesep, ND2files{i}];
counter =0; 

finfo = nd2finfo(filepath);
num_sites = finfo.img_seq_count/frames;  

for timept = 1:frames
for row = i
    for col = 1
        for site =1:num_sites
       
            
            shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
            
       if size (channels, 1) == 1
             [CFP]= nd2read_hayer_1chan(filepath,finfo, counter,1); 
             imwrite(CFP,[rawdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
        end 
            
       if size(channels, 1) ==2
         [FRET, CFP]= nd2read_hayer_2chan(filepath,finfo, counter,1); 
         imwrite(FRET,[rawdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
         imwrite(CFP,[rawdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
       end 
       
       if size (channels,1) == 3
              [mRuby, cyto,CAAX]= nd2read_hayer(filepath,finfo, counter,1); 
        imwrite(FRET,[rawdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(CFP,[rawdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[rawdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
       end 
            
       if size (channels,1) == 4
              [hoechst, pERM, mRuby, cyto]= nd2read_hayer(filepath,finfo, counter,1); 
           imwrite(hoechst,[rawdir,filesep,shot,'_hoechst_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(pERM,[rawdir,filesep,shot,'_pERM_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(mRuby,[rawdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
        imwrite(cyto,[rawdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
       end 
            
      if size (channels,1) == 6
             [mRuby, blank, cyto, blank2, FRET, CFP ] = nd2read_hayer_6chan(filepath,finfo,counter,1);
     imwrite(mRuby,[rawdir,filesep,shot,'_mRuby_',num2str(timept),'.tif'],'TIFF','Compression','None');
     imwrite(cyto,[rawdir,filesep,shot,'_cyto_',num2str(timept),'.tif'],'TIFF','Compression','None');
     imwrite(FRET,[rawdir,filesep,shot,'_FRET_',num2str(timept),'.tif'],'TIFF','Compression','None');
       imwrite(CFP,[rawdir,filesep,shot,'_CFP_',num2str(timept),'.tif'],'TIFF','Compression','None');
      end 
      %
      
      
    
      
      counter = counter+1;
      
            end   
        end 
    end 
end 
end 
 


%% tiff stacks from raw tiff images  
if makeTiffStacks == 1
    
    
 for row=1:(size(ND2files,1))
    
    for col = 1
        for site =1:10  % specify number of sites needed here 
            
          
           
writeTiffStacks(root, channel, row, col, site, frames); 
        end 
    end 
end   
end 

 %% generate averaged background images 

for chan=1:numel(channel1')
    calculate_bg_img_rm_blobs(bgdir,channel1{chan});
    disp(channel1(chan));
    disp(num2str(chan));
end
%% CFP/FRET channel alignment 
% when performing FRET imaging, small, microscope depedent shifts between
% the CFP and FRET channels can create artifacts. This section makes
% superimposed images of each channel and then generates alignment
% parameters that are used to align each channel to each other. 
if size (channels, 1)> 1 
CFPStack=[];
FRETStack=[];

k=0;
for row=1:4 % specify number of rows 
    
   
    for col=1 %specify number of columns in your raw data 

        for site=1:32 % specify number of sites
            k=k+1;
            shot=[num2str(row),'_',num2str(col),'_',num2str(site)];
            CFP_temp=imread([rawdir,filesep,shot,'_CFP_4.tif']); %choose which frame to use for alignment
            FRET_temp=imread([rawdir,filesep,shot,'_FRET_4.tif']); % CFP and FRET frame must be the same here! 
         
            CFPStack(:,:,k)=CFP_temp;
            FRETStack(:,:,k)=FRET_temp;
        
            
        end
    end
end


CFP_AV=uint16(mean(CFPStack,3));
FRET_AV=uint16(mean(FRETStack,3));

imwrite(CFP_AV,[bgdir,filesep,'AVG_rawdata_mRuby.tif'],'TIFF','Compression','None');
imwrite(FRET_AV,[bgdir,filesep,'AVG_rawdata_cyto.tif'],'TIFF','Compression','None');

alignStack(:,:,2)=imread([bgdir,filesep,'AVG_rawdata_cyto.tif']);
alignStack(:,:,1)=imread([bgdir,filesep,'AVG_rawdata_mRuby.tif']);
[pX,pY,dxMat1,dyMat1]=dualviewComputeAlignmentFromGridImages(alignStack);
figure;
subplot(1,2,1); imagesc(dxMat1); colorbar
subplot(1,2,2); imagesc(dyMat1); colorbar

save([bgdir,filesep,'alignment parameters pX pY.mat'],'pX','pY','dxMat1','dyMat1');
end 
%% generate raw FRET Data for full scale images (optional)
% necessary to generate bleaching correction curves for your channels, so 
% recommended to run anyways
% generates cell masks, cell centroid coordinates, as well as background corrected FRET
% ratios and channel intensities 

load([bgdir,filesep,'alignment parameters pX pY.mat']);

threshold = 4 ; % used for histogram based thresholding for mask generation 
k=0;
for row=1:2  %specify number of rows, columns, and sites.   
    for col=1:4
        for site=1:10


            k=k+1;
            
            position{k}=[num2str(row),'_',num2str(col),'_',num2str(site)];
        end
    end
end 

 
% if the threshold chosen is not correct for effective thresholding, just
% stop running the section, adjust the threshold, and rerun. 

%choose which function is appropriate for your needs. 
for k=1:length(position)
  % 1 channel 
     %getFRETDataHCS_1chan(position{k},bgdir,rawdir,datadir,threshold); 
% 2 channel FRET 
    getFRETDataHCS_2ChanFRET(position{k},bgdir,rawdir,datadir,threshold); 
% 3 channels: FRET and 1 extra channel, here listed as mRuby
    %getFRETDataHCS_3chan(position{k},bgdir,rawdir,datadir,threshold)
% 4 channels: CFP, FRET, mRuby, cytoplasm marker   
    %getFRETDataHCS_4chan(position{k},bgdir,rawdir,datadir,threshold)
end
disp('done!');
%clc; clear; 

%% bleaching correction
%generates bleaching decay curves which are used to normalize subsequent 
%single cell data

%FRET bleaching curve 
bleach_raw_all=[];
for k=1:length(position)
    load([datadir,filesep,position{k},'_Bleach_raw.mat']);
    bleach_raw_all(:,k)=bleach_raw/nanmedian(bleach_raw);
end
bleach_mean=nanmedian(bleach_raw_all(1:75, :),2);
plot(1:length(bleach_mean),bleach_raw_all(1:75, :));hold on
plot(1:length(bleach_mean),bleach_mean,'linewidth',3);

title('Bleaching Curve'); 
xlabel('TimePoint');
ylabel('Mean Intensity');

timepts=1:length(bleach_mean);
fitpara=fit(timepts',bleach_mean,'exp2');
corr=feval(fitpara,1:length(bleach_mean));
plot(1:length(bleach_mean),corr,'g');
axis([0 75 0.7 1.5]);
save([datadir,filesep,'bleachingcurve.mat'],'fitpara');


%mRuby bleaching 
bleach_raw_all=[]; for k=1:length(position)
    load([datadir,filesep,position{k},'_Bleach_raw.mat']);
    bleach_raw_all(:,k)=bleach_raw_mRuby/nanmedian(bleach_raw_mRuby);
end
bleach_mean=nanmedian(bleach_raw_all (1:75, :),2);
plot(1:length(bleach_mean),bleach_raw_all(1:75, :) );hold on
plot(1:length(bleach_mean),bleach_mean,'linewidth',3);

title('Bleaching Curve'); xlabel('TimePoint'); ylabel('Mean Intensity');

timepts=1:length(bleach_mean);
fitpara_mRuby=fit(timepts',bleach_mean,'exp2');
corr=feval(fitpara_mRuby,1:length(bleach_mean));
plot(1:length(bleach_mean),corr,'g'); axis([0 75 0 3]);

 save([datadir,filesep,'bleachingcurve_mRuby.mat'],'fitpara_mRuby');
 
 pause; 

 %cytoplasm marker bleaching
 bleach_raw_all=[]; for k=1:length(position)
    load([datadir,filesep,position{k},'_Bleach_raw.mat']);
    bleach_raw_all(:,k)=bleach_raw_cyto/nanmedian(bleach_raw_cyto);
end
bleach_mean=nanmedian(bleach_raw_all (1:75, :),2);
plot(1:length(bleach_mean),bleach_raw_all(1:75, :) );hold on
plot(1:length(bleach_mean),bleach_mean,'linewidth',3);

title('Bleaching Curve'); xlabel('TimePoint'); ylabel('Mean Intensity');

timepts=1:length(bleach_mean);
fitpara_cyto=fit(timepts',bleach_mean,'exp2');
corr=feval(fitpara_cyto,1:length(bleach_mean));
plot(1:length(bleach_mean),corr,'g'); axis([0 75 0 3]);

 save([datadir,filesep,'bleachingcurve_cyto.mat'],'fitpara_cyto');

%% corrected FRET data for full fields of view (optional) 
for k=1:length(position)
   
   %for CFP/YFP FRET only 
        correctBleachingExp2_FRET(position{k},fitpara, datadir); 
   
    % for CFP/YFP FRET, and 1 non ratio channel
        % correctBleachingExp2(position{k},fitpara, datadir, fitpara_mRuby);

   
end 

%% move to next section
open("single_cell_processing.m"); 
