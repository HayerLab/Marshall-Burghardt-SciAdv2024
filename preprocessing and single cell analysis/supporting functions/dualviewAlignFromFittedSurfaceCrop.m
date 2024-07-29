function im3=dualviewAlignFromFittedSurfaceCrop(im1,pX,pY,binning)

% Added a common cropping paratmeter to dualviewAlignFromFittedSurface.m, 
% so that the function can be run for multiple alignments 
% The resulting images are cropped symmetrically relative to the original images by 7
% pixels, right, left, top, bottom, irrespective of the shift.
% 210404 AH

cropsize=7; % number of pixels by which images are symmetrically cropped on all sides relative to the original image

if nargin<4
    %binning=round(1040/size(im1,1));
    binning=round(1080/size(im1,1)); % changed to match sCMOS chip and 1920x1080 capture size
end

%These values will be used to make the coordinate grid for the image being
%aligned
xVals1=1:size(im1,2);
yVals1=1:size(im1,1);

%Account for binning when applying the fitted function
unbinnedSize=size(im1(:,:,1))*binning;
xVals2=1:unbinnedSize(2);
yVals2=1:unbinnedSize(1);
xMat=dualviewGetXMatGridFor2ndOrderSurface(xVals2,yVals2);

fitdX=reshape([ones(unbinnedSize(1)*unbinnedSize(2),1) xMat]*pX,unbinnedSize); fitdX=imresize(fitdX,1/binning)/binning;
fitdY=reshape([ones(unbinnedSize(1)*unbinnedSize(2),1) xMat]*pY,unbinnedSize); fitdY=imresize(fitdY,1/binning)/binning;

[xCoor1,yCoor1]=meshgrid(xVals1,yVals1);
xCoor2=xCoor1-fitdX;   %This should give coordinates for image2 in the frame of image1
yCoor2=yCoor1-fitdY;
xCoor1rev=xCoor1+fitdX;  %This should give coordinates for image1 in the frame of image2
yCoor1rev=yCoor1+fitdY;

%Find the range of pixels in image1 that will be in the output image
yRange=[max(1,max(min(yCoor2))) min(size(im1,1),min(max(yCoor2)))];
yRange(1)=ceil(yRange(1));
yRange(2)=floor(yRange(2));
xRange=[max(1,max(min(xCoor2,[],2))) min(size(im1,2),min(max(xCoor2,[],2)))];
xRange(1)=ceil(xRange(1));
xRange(2)=floor(xRange(2));

%Find the output coordinates for the mapped image2, but in the frame of image2
xOut=xCoor1rev(yRange(1):yRange(2),xRange(1):xRange(2));
yOut=yCoor1rev(yRange(1):yRange(2),xRange(1):xRange(2));

im2(:,:,1)=im1(yRange(1):yRange(2),xRange(1):xRange(2),1);
im2(:,:,2)=interp2(xCoor1,yCoor1,im1(:,:,2),xOut,yOut);

%Crop by 'cropsize' number of pixels, relative to the original image
sizex_orig=size(im1,1); % size of original image, y
sizey_orig=size(im1,2); % size of original image, x
sizex_corr=size(im2,1); % size of shifted image, y
sizey_corr=size(im2,2); % size of shifted image, x

xRangeCrop(1)=cropsize-xRange(1)+2;
xRangeCrop(2)=sizex_orig-2*cropsize+xRangeCrop(1)-1;

yRangeCrop(1)=cropsize-yRange(1)+2;
yRangeCrop(2)=sizey_orig-2*cropsize+yRangeCrop(1)-1;

im3=im2(yRangeCrop(1):yRangeCrop(2),xRangeCrop(1):xRangeCrop(2),:);  


