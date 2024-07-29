function [maskFinal, cellCoors] = getCellMaskCyto_2(image, minCellSize, frameNum, position, threshold)
% getMask returns a binary mask based on image and minCellSize
% Uses a histogram-based thresholding
% Stores centroid and area of detected objects

%built in plotting option that can be used to visualize
%foreground-background separation

% Arnold Hayer 150520
% Seph Marshall June 18th, 2020




if nargin<2 %stands for number of argument inputs
    minCellSize=1500; % adjustable 
end


imSmooth1=imfilter(image,fspecial('disk',2),'replicate');
image_lpf=imfilter(imSmooth1,fspecial('disk',30),'symmetric');
imSmooth2=imSmooth1-image_lpf*0.25; % arbitrary scaling factor - adjust as needed



% Histogram-based threshold determination
% find range of pixel values
PixMinMax=double([round(prctile(imSmooth2(:),1)), round(prctile(imSmooth2(:),99))]); % can also use zero for the bottom limit 
IntStep=ceil((PixMinMax(2)-PixMinMax(1))/300);  

%plot probability dist. of pixel intensities in the frame
[f,xi]=ksdensity(imSmooth2(:),-500:PixMinMax(1,2)+1000); % can increase the range as much or as little as you want 
 
 a = strcmpi(position, '1_1_1');

if frameNum ==1 && a ==1
 figure; plot(xi,f); hold on;  
  end 
[pks,locs]=findpeaks(f,xi);

%% do this if cell edge is messy, lots of trailing bits 
%pks_cell = pks; 
%locs_cell =locs; 

%% filter peaks greater than 0.005
log_pks=pks>0.001; 
pks=pks(log_pks); locs=locs(log_pks);

 if frameNum ==1 && a ==1
        plot((locs),pks,'k^','markerfacecolor',[1 0 0]);
 end 
x_bgMax=locs(1,1); % picks the first peak (x-value of first peak)
[~,ind]=find(f>(0.01*pks(1)),1); % returns the first value of f greater than 0.5% of its max
x_1pct=xi(ind); % returns the corresponding intensity value
bgWidth=x_bgMax-x_1pct; % estimates the width of the background peak

 if frameNum ==1&& a ==1
      xline(x_bgMax+bgWidth); % to see how effective the estimation above of bg width is 
 end 
%Determine threshold for distinction between foreground/background (most important part of the code here)


threshSeg=(x_bgMax+(threshold)*bgWidth);% number here adjustable: if having trouble with segmentation adjust based on fg/bg separation%xline(threshSeg,'--'); % again just for visualization 

if frameNum ==1 && a ==1
          xline(threshSeg, '--'); 
          pause; 
                 hold off;
              
 end 

%% messy cell edge part 2 
% uses average cell intensity instead of bg peak, and thresholds backwards
% from this 
%locs_cell_boolean=locs_cell>threshSeg;
%locs_cell=locs_cell .*(locs_cell_boolean); 
%pks_cell= pks_cell.*(locs_cell_boolean); 

%finds cell peak through increased AROC on pixel intensity graph
% for i =10:size(pks_cell,2)
%     
%   if (pks_cell(1,i)-pks_cell(1,i-9)) >= 0.001
%       break; 
%   end 
%     
% end 

%[value, index] = max(pks_cell); 


 
%plot(locs_cell(1,index), pks_cell(1,index), 'k^','markerfacecolor',[1 0 0]);


%threshSeg=(locs_cell(1,index)/1.8); % adjustable parameter
%xline(threshSeg,'r:')

%hold off; 
%% masks the picture, then removes small areas of
%background less than 200 pixels surrounded by cell (usually mistakes);

mask_init=imSmooth2>threshSeg;
%imagesc(mask_init); %for visualization if necessary 
mask=bwareaopen(mask_init,minCellSize); 
%imagesc(mask); %for visualization if necessary 



 background = ~mask; 
  background = bwareaopen(background, 400); % added this to take out little flickers inside cells 
 mask_filled1 = ~background; % added this to take out little flickers inside cells 
%imagesc(mask_filled1); %for visualization if necessary 

% this section is for if masking a monolayer with holes in between cells 

%       holes=mask_filled1 & ~mask;
%       bigholes =bwareaopen(holes,50);
%     smallholes= holes & ~bigholes;
% %     
    % maskFinal = mask | smallholes;

mem_mask2 = imopen(mask_filled1,strel('disk',2)); %get rid of tiny fibers
maskFinal = mask_filled1 & mem_mask2; 

maskFinal= bwareaopen(maskFinal,minCellSize);
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