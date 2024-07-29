function [maskFinal, cellCoors] = getCellMaskCyto_2_stacked(image, minCellSize, frameNum,s, threshold)
% getMask returns a binary mask based on image and minCellSize
% Uses a histogram-based thresholding
% Stores centroid and area of detected objects

%built in plotting option that can be used to visualize
%foreground-background separation

% Arnold Hayer 150520
% Seph Marshall June 18th, 2020




if nargin<2 %stands for number of argument inputs
    minCellSize=1000; % adjustable 
end


imSmooth1=imfilter(image,fspecial('disk',1),'replicate');
image_lpf=imfilter(imSmooth1,fspecial('disk',30),'symmetric');
imSmooth2=imSmooth1-image_lpf*0.25; % arbitrary scaling factor - adjust as needed



% Histogram-based threshold determination
% find range of pixel values
PixMinMax=double([round(prctile(imSmooth2(:),1)), round(prctile(imSmooth2(:),99))]); % can also use zero for the bottom limit 
IntStep=ceil((PixMinMax(2)-PixMinMax(1))/300);  

%plot probability dist. of pixel intensities in the frame
[f,xi]=ksdensity(imSmooth2(:),-500:PixMinMax(1,2)+1000); % can increase the range as much or as little as you want 
 
 

if frameNum ==1 && s ==1
 figure; plot(xi,f); hold on;  
end 
[pks,locs]=findpeaks(f,xi);

%% do this if cell edge is messy, lots of trailing bits 
%pks_cell = pks; 
%locs_cell =locs; 

%% filter peaks greater than 0.005
log_pks=pks>0.001; 
pks=pks(log_pks); locs=locs(log_pks);

 if frameNum ==1 && s ==1
        plot((locs),pks,'k^','markerfacecolor',[1 0 0]);
 end 
 

x_bgMax=locs(1,1); % picks the first peak (x-value of first peak)
[~,ind]=find(f>(0.01*pks(1)),1); % returns the first value of f greater than 0.5% of its max
x_1pct=xi(ind); % returns the corresponding intensity value
bgWidth=x_bgMax-x_1pct; % estimates the width of the background peak

 if frameNum ==1  && s ==1
      xline(x_bgMax+bgWidth); % to see how effective the estimation above of bg width is 
 end 
%Determine threshold for distinction between foreground/background (most important part of the code here)


threshSeg=(x_bgMax+(threshold)*bgWidth);% number here adjustable: if having trouble with segmentation adjust based on fg/bg separation%xline(threshSeg,'--'); % again just for visualization 

if frameNum ==1  && s ==1
          xline(threshSeg, '--'); 
          pause; 
                 hold off;
              
 end 




%% masks the picture

mask_init=imSmooth2>threshSeg;
%imagesc(mask_init); %optional visualization
mask=bwareaopen(mask_init,minCellSize); 
%imagesc(mask); %optional visualization

mask = bwareafilt(mask,1); % chooses the largest mask in the image - excludes small mask parts from other cells in periphery
%imagesc(mask); %optional visualization



background = ~mask; 
background = bwareaopen(background, 400); % added this to take out little flickers inside cells 
mask_filled1 = ~background; % added this to take out little flickers inside cells 
%imagesc(mask_filled1); 

% this section is for if masking a monolayer with holes in between cells 

%     holes=mask_filled1 & ~mask;
%     bigholes =bwareaopen(holes,1000);
%     smallholes= holes & ~bigholes;
%     
%     maskFinal = mask | smallholes;

mem_mask2 = imopen(mask_filled1,strel('disk',2)); %get rid of unwanted tiny fibers
maskFinal = mask_filled1 & mem_mask2; 

maskFinal= bwareaopen(maskFinal,minCellSize); %thisMask
%imagesc(maskFinal);

% get centroid coordinates
objects=regionprops(maskFinal,'BoundingBox','Area','Centroid','PixelIdxList','PixelList');
objects=objects(arrayfun(@(x) x.Area>minCellSize,objects));
%objects=objects(arrayfun(@(x) x.Area<maxCellSize,objects));

cellCoors(:,1)=vect(arrayfun(@(x) x.Centroid(1),objects));
cellCoors(:,2)=vect(arrayfun(@(x) x.Centroid(2),objects));
cellCoors(:,3)=vect(arrayfun(@(x) x.Area,objects));
%save('maskFinal');
end