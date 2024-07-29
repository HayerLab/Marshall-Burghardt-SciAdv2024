function [ ActinMask ] = getActinMaskIF(im_in,lowpix)
% getActinMaskIF returns a binary mask based on images of subconfluent to near-confluent, 
% phalloidin-stained cells
% Uses the 0.02% of pixels to estimate the background, mask is based on
% a threshold equivalent to 1.3* median of background pixel intensity.  
% User-defined percentile of pixels considered for mask  (lowpix) , default is 0.02%
% Arnold Hayer 230521

if nargin<2
    percentile=0.02;
end 

im_in_filt=imfilter(im_in,fspecial('disk',5),'replicate');
allpix=im_in_filt(:);
lowintpixels=allpix(allpix<prctile(allpix,lowpix)); %identifies 0.02% lowest pixles (1062 pixels for 2304x2304)
thrsh=median(lowintpixels)*1.3;
ActinMask=im_in_filt>thrsh;

%subplot(1,2,1); imagesc(im_in_filt);
%subplot(1,2,2); imagesc(ActinMask);
end


