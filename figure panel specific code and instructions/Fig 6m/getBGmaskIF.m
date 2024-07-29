function [ BGmask ] = getBGmaskIF(im_in,lowpix)
% getBGmaskIF a binary mask based on image for background subtraction
% Background pixels are assumed to be the ones with the lowest intensities, 
% User-defined percentile of pixels considered for mask  (lowpix) , default is 0.5%
% Arnold Hayer 230521

if nargin<2
    lowpix=0.5;
end 

im_in_filt=imfilter(im_in,fspecial('disk',3),'replicate');
lowintensitypixels=prctile(im_in_filt(:),lowpix);
thrsh=lowintensitypixels;
BGmask=im_in_filt>thrsh;
end
