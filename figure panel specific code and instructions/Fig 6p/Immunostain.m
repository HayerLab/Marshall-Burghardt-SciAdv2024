%% Immunostain
% To analyze IF signal in immunostained HUVEC
% 
% pERM quantification
% - Based on cytoplasmic rings (ringsig3_med), or expanded nuclear masks (expnucsig3_med) surrounding Hoechst-defined nuclei.
% - Paramters set to work with images acquired using a 20x lens, no binning.
% - Background levels are based on bg images and a BG mask 
% - Output is in expnucsig3_med or 
%
% pMLC
% - Based on a phalloidin mask, then quanitificaiton of pMLC signal within
% the mask, (sig3FOV_mean). 
%
% Input: 
%
% IF images are labeled row_col_site_channel.tif (shot_channel.tif)
% 
% Output: 
% 
% Per image, one output file is generated (IFData_shot.mat]), in datadir, containing
% an array IFData: 
    
    % Column 1: x-coordinate of cell nucleus.
    % Column 2: y-coordinate of cell nucleus.
    % Column 3: nuclear area (based in signal 1)
    % Column 4: mean signal 2 intensity within cytoplasmic ring (ringwidth) of a given cell.
    % Column 5: mean signal 3 intensity within cytoplasmic ring (ringwidth) of a given cell.
    % Column 6: mean signal 2 intensity within expanded nucleus (ringwidth)
    % Column 7: mean signal 3 intensity within expanded nucleus (ringwidth)
    % Column 8: mean signal 3 intensity within mask based on signal 2

% Arnold Hayer, 21 May 2023

clear;clc;close all;

root='..\230714_TI2E';
rawdir=[root,filesep,'Raw'];
bgpath='..\background';
datadir=[root,filesep,'data'];
if ~exist(datadir)
    mkdir(datadir);
end

%% Background images and parameters
bgim1=single(imread([bgpath,filesep,'bg_395.tif']));
bgim2=single(imread([bgpath,filesep,'bg_470.tif']));
bgim3=single(imread([bgpath,filesep,'bg_555.tif']));

name1='395 nm'; %nuc
name2='470 nm'; %phalloidin-488
name3='555 nm'; %anti-pERM-568, anti-pMLC

ringwidth=7; % width of the nuclear ring/nuclear expansion in pixels

moviebin=1;

nucedgename='nucedge';

if moviebin==1
    nucr=24; 
    debrisarea=800; 
    
elseif moviebin==2
    nucr=6;
    debrisarea=50; 
end
boulderarea=20*debrisarea; 
blobthreshold=-0.06; 
timetotal=tic;

%% Loop through image files
for rows=2:7
    for cols=2:11
        for sites=1:16
            shot=[num2str(rows),'_',num2str(cols),'_',num2str(sites)];
            disp(shot);
            IFdata=[];
            try % try-catch is used here to avoid code from stopping in case wells/images were removed for analysis (remove for debugging)
                    im1raw=single(imread([rawdir,filesep,shot,'_',name1,'_1.tiff'])); 
                    im2raw=single(imread([rawdir,filesep,shot,'_',name2,'_1.tiff']));
                    im3raw=single(imread([rawdir,filesep,shot,'_',name3,'_1.tiff']));
                    % segment nuclei 
                    nuc_mask=blobdetector(log(im1raw),nucr,blobthreshold,debrisarea);
                    foreground=nuc_mask;
                    nuc_mask=segmentdeflections(nuc_mask,nucr,0.5,debrisarea); %comment out for sparse YT analysis
                    nuc_mask=excludelargeandwarped(nuc_mask,boulderarea);
                    % remove border objects 
                    [height,width]=size(im1raw);
                    nuc_mask([1 height],1:width)=1; nuc_mask(1:height,[1 width])=1;
                    border=bwareaopen(nuc_mask,height*2+width*2-4);
                    nuc_mask=logical(nuc_mask-border);
                    % Subtract background for name 1, name2, and name3
                    im1filt=imfilter(im1raw,fspecial('disk',3),'symmetric');
                    mask1=getBGmaskIF(im1filt,5);% 0.5 for close fo confluent cultures 5 for sparse
                    im1bgsub=subBG(im1filt,mask1,bgim1);

                    im2filt=imfilter(im2raw,fspecial('disk',3),'symmetric');
                    mask2=getBGmaskIF(im2filt,0.5);%0.5 for close fo confluent cultures 5 for sparse
                    im2bgsub=subBG(im2filt,mask2,bgim2);

                    im3filt=imfilter(im3raw,fspecial('disk',3),'symmetric');
                    mask3=getBGmaskIF(im3filt,0.5); %0.5 for close fo confluent cultures 5 for sparse
                    im3bgsub=subBG(im3raw,mask3,bgim3);

                    nuc_info=struct2cell(regionprops(nuc_mask,'Area','Centroid')');
                    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
                    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';

                % extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                nuc_label=bwlabel(nuc_mask);
                numcells=numel(nuc_area);
                nuc_info=regionprops(nuc_label,'PixelIdxList');
                %ring_label=getcytoring(nuc_label);
                ring_label=getcytoring_3(nuc_label,ringwidth,im3bgsub); %10xB1:4 20xB1:4
                ring_info=regionprops(ring_label,'PixelIdxList');
                
                exp_nuc_label=getexpandnuc(nuc_label, ringwidth,im3bgsub);
                exp_nuc_info=regionprops(exp_nuc_label,'PixelIdxList');
                
                nanvec=ones(numcells,1)*NaN;
                sig2ring=nanvec;
                sig3ring=nanvec;
                
                sig2expnuc2=nanvec;
                sig3expnuc3=nanvec;

                for cc=1:numcells % loop through objects
                        %sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
                        %sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
                        %sig3(cc)=median(real3(nuc_info(cc).PixelIdxList));
                        %sig2ring_75th(cc)=prctile(ringall,75); %previously mean
                        %sig1(cc)=mean(im1bgsub(nuc_info(cc).PixelIdxList));

                        ringall2=im2bgsub(ring_info(cc).PixelIdxList); 
                        ringall2(ringall2>prctile(ringall2,95))=[]; %removes 5% brightest pixels;
                        sig2ring(cc)=mean(ringall2);
                        
                        expnucall2=im2bgsub(exp_nuc_info(cc).PixelIdxList);
                        expnucall2(expnucall2>prctile(expnucall2,95))=[];
                        sig2expnuc2(cc)=mean(expnucall2);


                        ringall3=im3bgsub(ring_info(cc).PixelIdxList); 
                        ringall3(ringall3>prctile(ringall3,95))=[]; %removes 5% brightest pixels;
                        sig3ring(cc)=mean(ringall3);
                        
                        expnucall3=im3bgsub(exp_nuc_info(cc).PixelIdxList);
                        expnucall3(expnucall3>prctile(expnucall3,95))=[];
                        sig3expnuc3(cc)=mean(expnucall3);
                end

                % pMLC
                sig3FOV=[];
                sig2mask=getActinMaskIF(im2raw,0.02); %mask based on sig2 (actin-phalloidin)
                imagesc(sig2mask);
                sig3FOV=median(vect(im3bgsub(sig2mask))); %sig3 (pMLC) in sig2mask

                %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1];
                %imID=ones(numcells,1).*frame; %to track back in the array individual images
                IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig2ring,sig3ring,sig2expnuc2,sig3expnuc3];
                %IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3];
                %IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3,sig4];
                %IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig2ring_75th,sig2ring_50th,sig2ring_nobg,sig3,sig4];
                %IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1,sig2,sig3,sig4,sig2_ring,sig3_ring];
                save([datadir,filesep,'IFdata_',shot,'.mat'],'IFdata','sig3FOV');
            end %try
            end %sites
    end %frames
end %rows

disp('done!');

%% Compute averaged data 
% averaging per well. 

clear;clc;close all;
root='..\240417_TI2E';
datadir=[root,filesep,'data'];

for rows=2:7
    for cols=2:11
        ringsig3=[];
        expnucsig3=[];
        sig3FOV_sites=[];
        for sites=1:16
            shot=[num2str(rows),'_',num2str(cols),'_',num2str(sites)];
            disp(shot);
            try
               load([datadir,filesep,'IFdata_',shot,'.mat']);
               

               ringsig3=[ringsig3; IFdata(:,5)];
               expnucsig3=[expnucsig3;IFdata(:,7)];
               sig3FOV_sites=[sig3FOV_sites;sig3FOV];
            end
        end
      ringsig3_plate{rows,cols}=ringsig3;
      expnucsig3_plate{rows,cols}=expnucsig3;
      ringsig3_med(rows,cols)=median(ringsig3);
      expnucsig3_med(rows,cols)=median(expnucsig3);
      sig3FOV_mean(rows,cols)=mean(sig3FOV_sites);
    end
end
disp('done!');

