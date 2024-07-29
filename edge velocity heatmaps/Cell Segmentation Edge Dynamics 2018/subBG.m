function [ im_out ] = subBG(im_in,mask,bg_image )
% subtracts the median background value as determined from the backgroud
% mask
% 
if nargin<3
    bg=nanmedian(vect(im_in(~mask)));
    im_out=im_in-bg;
else 
    bg_level=nanmedian(vect(im_in(~mask))); % median intensity of background within image
    bg_bg=nanmedian(vect(bg_image(~mask))); % median intensity within mask in background image
    bg=bg_image.*(bg_level/bg_bg); % normalization background image to actual background levels
    im_out=im_in-bg;
end
end

