re%% single cell cell preprocessing before edge tracking
% manually crops cells of interest from tiff_stacks 
%generates cell masks, and FRET/fluorescent channel data 

%% initiatialization
clear; clc; 

 root = 'F:\example_dataset_240522';
 
 tiffDir = ([root, filesep,'tiff_stacks']); 
    
     
cropdir=[root,filesep,'cropped'];
if ~exist(cropdir)
         mkdir(cropdir); 
end  
bgpath=[root,filesep,'background'];
   
    %% crop cell files and accompanying background 

% run this section as many times as needed, just change the cell number 
% which specifies which folder to save under and filekey which specifies 
% which orginal tiff stack to draw from 
   clc; 
cell=1; 
filekey = '1_1_1'; 


 cellDir = ([cropdir,filesep, num2str(cell)]); 
     if ~exist(cellDir)
         mkdir(cellDir); 
     end   
     
% choose which channels you are working with here. Currently configured for
% FRET, CFP, mRuby, and cytoplasm marker. Any can be commented out or
% removed if not needed

bgcyto = ([bgpath,filesep, 'AVG_bgcyto.tif']); 
bgCFP = ([bgpath, filesep, 'AVG-BG-CFP.tif']); 
bgFRET = ([bgpath,filesep,  'AVG_bgFRET.tif']);
bgmRuby=([bgpath, filesep,  'AVG_bgmRuby.tif']);
          

 CFP= ([tiffDir,filesep, strcat(filekey,'_CFP_stacked.tif')]); 
 FRET= ([tiffDir,filesep, strcat(filekey,'_FRET_stacked.tif')]);
 mRuby= ([tiffDir,filesep, strcat(filekey,'_mRuby_stacked.tif')]); 
 cyto= ([tiffDir,filesep, strcat(filekey,'_cyto_stacked.tif')]);  
  




 bg_CFP_image = double(readTIFFstack(bgCAAX)); 
 bg_FRET_image = double(readTIFFstack(bgFRET)); 
bg_mRuby_image = double(readTIFFstack(bgmRuby)); 
bg_cyto_image = double(readTIFFstack(bgcyto)); 
     
            cropSite = 0;


% here using mRuby channel to visualize the data and choose Crop sites. 
% all other channels and background images are cropped to the same
% dimensions. if you wish to use a different channel for cropping, simply
% switch the channel names 

mRuby_stack=double(readTIFFstack(mRuby));
         
        mRuby_stack_1 = mRuby_stack(:,:,1); 
               
         fg = figure;
              
         axis image;  
             
            
      [stackmRuby, cropArea] = serimcropold(mRuby_stack,mean(mRuby_stack,3));          
      Stack2TIFF(stackmRuby, [cellDir, filesep,'mRuby.tif']);
                    
       stackCFP = readFileToStack(CFP); 
       stackCFP = imcrop3(stackCFP, [cropArea(1), cropArea(2), 1,...
       cropArea(3), cropArea(4), size(stackCAAX,3)-1]);                    
       Stack2TIFF(stackCFP, [cellDir, filesep, 'CFP.tif']);
       
       stackFRET = readFileToStack(pERM); 
       stackFRET = imcrop3(stackFRET, [cropArea(1), cropArea(2), 1,...
       cropArea(3), cropArea(4), size(stackFRET,3)-1]);                    
       Stack2TIFF(stackFRET, [cellDir, filesep, 'FRET.tif']);

        stackcyto = readFileToStack(cyto); 
       stackcyto = imcrop3(stackcyto, [cropArea(1), cropArea(2), 1,...
       cropArea(3), cropArea(4), size(stackcyto,3)-1]);                    
       Stack2TIFF(stackcyto, [cellDir, filesep, 'cyto.tif']);
       
            
       

    
      FRET_bg_crop = imcrop(uint16(bg_FRET_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]);
      imwrite(FRET_bg_crop,[cellDir, filesep, 'FRET_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
                    
      CFP_bg_crop = imcrop(uint16(bg_CFP_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
     imwrite(CFP_bg_crop,[cellDir, filesep, 'CFP_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
          
      mRuby_bg_crop = imcrop(uint16(bg_mRuby_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]);              
      imwrite(mRuby_bg_crop,[cellDir, filesep, 'mRuby_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
            
            
      cyto_bg_crop = imcrop(uint16(bg_cyto_image), [cropArea(1), cropArea(2), cropArea(3), cropArea(4)]); 
      imwrite(cyto_bg_crop,[cellDir, filesep, 'cyto_bg.tif'] , "WriteMode", "overwrite", "Compression", "none");
                    
                   
            clear stack;
  close all; 
  
  disp('ready for next crop'); 
    
    

%% FRET data 
%loading of alignment parameters and bleaching correction data 

clc;


bleachdir=([root,filesep,'data']);
bgdir=[root,filesep,'background'];
  load([bleachdir,filesep,'bleachingcurve.mat']);
  load([bleachdir,filesep,'bleachingcurve_mRuby.mat']);
 load([bgpath, filesep, 'alignment parameters pX pY.mat']); 
 load([bleachdir,filesep,'bleachingcurve_cyto.mat']);


%% Parallel loop
% specify k in the for loop as number of cells you have cropped
threshold = 3; 
for k=1:29
  
   rawdir=[cropdir,filesep, strcat( num2str(k))]; 
    load([rawdir,filesep,'alignment parameters pX pY.mat']);
    
   datadir=[rawdir,filesep,'output'];

 if ~exist(datadir)
    mkdir(datadir);
 end
 
 %choose which one you need 
 
 % this one has FRET/CFP configured, commented out are options for a 3rd
 % and 4th channel if you want 
getFRETDataHCS_stacked_FRET(k,rawdir,datadir,threshold, pX, pY); %FRET, CFP
 %getFRETDataHCS_stacked_3chan(k,rawdir,datadir,3.5, pX, pY); % FRET, CFP, mRuby
 %getFRETDataHCS_stacked_4_chan(k,rawdir,datadir,1.7, pX, pY); % FRET, CFP,
 %mRuby, cyto used in RhoB sensor ezrin ratio 
%getFRETDataHCS_stacked_ezrin_cyto_caax(k,rawdir,datadir,3) % for ezrin
%ratio calculations, also has option for extra channel, in this case a membrane marker 


% choose which one you need

%correctBleachingExp2_stacked(fitpara, datadir, fitpara_mRuby); %  % does FRET, and mRuby
correctBleachingExp2_FRET_stacked(fitpara, datadir); %only does FRET
%correctBleachingExp2_cyto_ratio_stacked(datadir, fitpara_mRuby, fitpara_cyto); % for ezrin ratio calculations 
 
    
end
disp('done!');


