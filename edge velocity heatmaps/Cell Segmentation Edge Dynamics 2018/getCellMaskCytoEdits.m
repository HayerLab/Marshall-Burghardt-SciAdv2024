function [maskFinal cellCoors] = getCellMaskCytoEdits(image,  minCellSize)
% getMask returns a binary mask based on image and minCellSize
% Uses a histogram-based thresholding
% Stores centroid and area of detected objects
% Arnold Hayer 150520



% if nargin<3
%     maxCellSize=16000;
% end
if nargin<2
    minCellSize=2000;
end

% Membrane-stain: Subtract scaled blurred image to sharpen outline

imSmooth1=imfilter(image,fspecial('disk',1),'replicate');
image_lpf=imfilter(imSmooth1,fspecial('disk',30),'symmetric');
imSmooth2=imSmooth1-image_lpf*0.25; % arbitrary scaling factor

% Cytoplasmic stain: Subtract scaled blurred image to sharpen outline
% imSmooth2=imfilter(image,fspecial('disk',5),'replicate');

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

[m,n] =size(locs);

%plot(xi(locs),pks,'k^','markerfacecolor',[1 0 0]);hold off;

x_bgMax=xi(locs(1)); % picks the first peak (x-value of first peak)

[~,ind]=find(f>(0.01*pks(1)),1); % returns the first value of f greater than 0.5% of its max
x_1pct=xi(ind); % returns the corresponding intensity value
bgWidth=x_bgMax-x_1pct; % estimates the width of the background peak

if (n<2)
    x_cellMax=0;
    x_2pct=0; 
    cellWidth = 0;
    
    threshSeg=(x_bgMax+bgWidth);
else
    x_cellMax=xi(locs(2)); % here didnt attempt to find the width here, just first peak of cell fluorescence
  
    threshSeg = ((x_cellMax)+(x_bgMax+bgWidth))/2;
end 

 

%figure, showImagesWithLinkedAxes(imSmooth1,imerode(imSmooth2>threshSeg,strel('disk',1)));

%seph - removing the 2 here seems to make the edge outlining sharper 

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
