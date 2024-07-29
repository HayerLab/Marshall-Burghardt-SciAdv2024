function [protrusions,retractions,image]=getEdgeVelStats(inputmap,protThresh,retrThresh,sizeThresh)
% Input: 
%    inputmap - spatial/temporal edgevelocity map
%    retrThres - edge velocity threshold for protrusions
%    sizeThres - size threshold for protrusions
% Output:
%     Image
%     protrusions, retractions: Structure arrays with the following parameers: 
%     - Area
%     - Location (centroid)
%     - Image of protrusion (logical, with smallest possible outline)
%     - PixelValues
%     - x - duration
%     - y - size in number of windows
%     - median velocity. 
%     - Image with protrusions/retractions outlined. 

% Arnold Hayer, July 26, 2017

% retrThresh=-5;
% protThresh=5;

%% Filter out small junk
LOprot=inputmap>protThresh;LOprot=bwareaopen(LOprot,sizeThresh);protLabel=bwlabel(LOprot);
LOretr=inputmap<retrThresh;LOretr=bwareaopen(LOretr,sizeThresh);retrLabel=bwlabel(LOretr);
 %subplot(1,2,1);imagesc(LOprot);
 %subplot(1,2,2);imagesc(LOretr);
im1=ratio2RGB(inputmap,[-12 12]);
image=DrawMaskOutline(im1,LOprot+LOretr,[1 0 0]);
 imshow(image)

% protrusions
protrusions=regionprops(LOprot,inputmap,'Area','Centroid','Image','PixelValues');

for k=1:size(protrusions,1)
    imagesc(protrusions(k).Image);
    protrusions(k).xmax=max(sum(protrusions(k).Image,1));
    protrusions(k).ymax=max(sum(protrusions(k).Image,2));
    protrusions(k).mspeed=median(protrusions(k).PixelValues);
    pause;
end

% retractions
retractions=regionprops(LOretr,inputmap,'Area','Centroid','Image','PixelValues');

for k=1:size(retractions,1)
    imagesc(retractions(k).Image);
    retractions(k).xmax=max(sum(retractions(k).Image,1));
    retractions(k).ymax=max(sum(retractions(k).Image,2));
    retractions(k).mspeed=median(retractions(k).PixelValues);
    pause; 
end

end

