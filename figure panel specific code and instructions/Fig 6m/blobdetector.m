function nuc_mask = blobdetector_3(nuc_raw,nucr,threshold,debrisarea)
nuc_raw_sharp=imsharpenSean(nuc_raw,2,nucr*0.25); %added this on 2014/01/02
%nuc_raw_sharp=nuc_raw;
sigma=0.75*nucr/sqrt(2); %default 0.75
h=sigma^2*fspecial('log',[nucr*2 nucr*2],sigma); %laplacian of gaussian default window [nucr*2 nucr*2]
nuc_log=imfilter(nuc_raw_sharp,h,'symmetric');
nuc_mask=nuc_log<threshold; %higher picks up debris, lower misses parts of nuclei
nuc_mask=imfill(nuc_mask,'holes');
nuc_mask=imopen(nuc_mask,strel('disk',2,0)); %bin1:4 bin2:2. AH: was 2, changed to 1
nuc_mask=~bwmorph(~nuc_mask,'diag');
nuc_mask=~bwmorph(~nuc_mask,'bridge');
nuc_mask=bwareaopen(nuc_mask,debrisarea); %bin1:0.75 bin2:0.5
end

%%% debugging %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
extractmask=bwmorph(nuc_mask,'remove');
tempframe=imadjust(mat2gray(nuc_raw));
tempframe(:,:,2)=extractmask;
tempframe(:,:,3)=0;
figure,imshow(tempframe);
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function imSharp = imsharpenSean(im,magnitude,gaussianRadius,threshold)
%IMSHARPENSEAN implements image sharpening
%   Detailed explanation goes here

imBlurred=imfilter(im,fspecial('gaussian',max(3,round(2.5*gaussianRadius)),gaussianRadius),'symmetric');

sharpMask=im-imBlurred; %imagesc(sharpMask); colorbar;
imSharp=im+magnitude*sharpMask;
%imSharp=im-magnitude*imBlurred;
%showImagesWithLinkedAxes({im,imSharp});

end