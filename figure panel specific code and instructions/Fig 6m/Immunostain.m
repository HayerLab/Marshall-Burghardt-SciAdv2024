%% Immunostain
% To analyze IF signal in immunostained HUVEC
% 
% Immunosignal quantification
% - Based on cytoplasmic rings (ringsig2_med), or expanded nuclear masks (expnucsig2_med) surrounding Hoechst-defined nuclei.
% - Paramters set to work with images acquired using a 20x lens, no binning.
% - Background levels are based on bg images and a BG mask 
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
    % Column 5: mean signal 2 intensity within cytoplasmic ring (ringwidth) of a given cell.
    

% Arnold Hayer, 21 May 2023

clear;clc;close all;

root='..\240512_TI2E';
rawdir=[root,filesep,'Raw'];
bgpath='..\background';
datadir=[root,filesep,'data'];
if ~exist(datadir)
    mkdir(datadir);
end

%% Background images and parameters
bgim1=single(imread([bgpath,filesep,'bg_395.tif']));
bgim2=single(imread([bgpath,filesep,'bg_470.tif']));


name1='395 nm'; %nuc
name2='470 nm'; %phalloidin-488


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
                    if size(im1raw,1)==1152 % scale images 2x2 in case they were acquired at 2x2 binning
                        im1raw=imresize(im1raw,2)./4; % divided by 4 to normalize to intensity of no bin images
                        im2raw=imresize(im2raw,2)./4; % dividd by 4 to normalized to intensity of no bin images
                        
                    end
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
                    % Subtract background for name 1, name
                    im1filt=imfilter(im1raw,fspecial('disk',3),'symmetric');
                    mask1=getBGmaskIF(im1filt,5);% 0.5 for close fo confluent cultures 5 for sparse
                    im1bgsub=subBG(im1filt,mask1,bgim1);

                    im2filt=imfilter(im2raw,fspecial('disk',3),'symmetric');
                    mask2=getBGmaskIF(im2filt,0.5);%0.5 for close fo confluent cultures 5 for sparse
                    im2bgsub=subBG(im2filt,mask2,bgim2);

                    nuc_info=struct2cell(regionprops(nuc_mask,'Area','Centroid')');
                    nuc_area=squeeze(cell2mat(nuc_info(1,1,:)));
                    nuc_center=squeeze(cell2mat(nuc_info(2,1,:)))';

                % extract features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                nuc_label=bwlabel(nuc_mask);
                numcells=numel(nuc_area);
                nuc_info=regionprops(nuc_label,'PixelIdxList');
                %ring_label=getcytoring(nuc_label);
                ring_label=getcytoring_3(nuc_label,ringwidth,im2bgsub); %10xB1:4 20xB1:4
                ring_info=regionprops(ring_label,'PixelIdxList');
                
                exp_nuc_label=getexpandnuc(nuc_label, ringwidth,im2bgsub);
                exp_nuc_info=regionprops(exp_nuc_label,'PixelIdxList');
                
                nanvec=ones(numcells,1)*NaN;
                sig2ring=nanvec;
                sig3ring=nanvec;
                
                sig2expnuc2=nanvec;
                %sig3expnuc3=nanvec;

                for cc=1:numcells % loop through objects
                        %sig1(cc)=mean(real1(nuc_info(cc).PixelIdxList));
                        %sig2(cc)=median(real2(nuc_info(cc).PixelIdxList));
                        %sig2ring_75th(cc)=prctile(ringall,75); %previously mean
                        %sig1(cc)=mean(im1bgsub(nuc_info(cc).PixelIdxList));

                        ringall2=im2bgsub(ring_info(cc).PixelIdxList); 
                        ringall2(ringall2>prctile(ringall2,95))=[]; %removes 5% brightest pixels;
                        sig2ring(cc)=mean(ringall2);
                        
                        expnucall2=im2bgsub(exp_nuc_info(cc).PixelIdxList);
                        expnucall2(expnucall2>prctile(expnucall2,95))=[];
                        sig2expnuc2(cc)=mean(expnucall2);
                       
                end

                %%% store data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig1];
                %imID=ones(numcells,1).*frame; %to track back in the array individual images
                IFdata=[nuc_center(:,1),nuc_center(:,2),nuc_area,sig2ring,sig2expnuc2];
                save([datadir,filesep,'IFdata_',shot,'.mat'],'IFdata');
            end %try
            end %sites
    end %frames
end %rows

disp('done!');

%% Compute averaged data 
% averaging per well. 

clear;clc;close all;
root='..\240512_TI2E';
datadir=[root,filesep,'data'];

for rows=2:7
    for cols=2:11
        ringsig2=[];
        expnucsig2=[];
        for sites=1:16
            shot=[num2str(rows),'_',num2str(cols),'_',num2str(sites)];
            disp(shot);
            try
               load([datadir,filesep,'IFdata_',shot,'.mat']);
               ringsig2=[ringsig2; IFdata(:,4)];
               expnucsig2=[expnucsig2;IFdata(:,5)];

            end
        end
     
      ringsig2_plate{rows,cols}=ringsig2;
      expnucsig2_plate{rows,cols}=expnucsig2;
      ringsig2_med(rows,cols)=median(ringsig2);
      expnucsig2_med(rows,cols)=median(expnucsig2);
    end
end
disp('done!');

