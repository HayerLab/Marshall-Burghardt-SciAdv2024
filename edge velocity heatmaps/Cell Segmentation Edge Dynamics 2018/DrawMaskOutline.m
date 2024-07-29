function [ outlined_mask ] = DrawMaskOutline( ImIn,mask,linecolor )
% DrawMaskOutline draws a green outline around masked objects in a greyscale or RGB image image
% based on a mask onto an image
% Arnold, 150521
if nargin<3
    linecolor=[0,1,0];
end

BWoutline = bwperim(mask);
if size(ImIn,3)==1
    ImIn=imadjust(mat2gray(ImIn));
    redChan=ImIn; redChan(BWoutline)=linecolor(1);
    greenChan=ImIn; greenChan(BWoutline)=linecolor(2);
    blueChan=ImIn; blueChan(BWoutline)=linecolor(3);
    outlined_mask(:,:,1)=redChan;
    outlined_mask(:,:,2)=greenChan;
    outlined_mask(:,:,3)=blueChan;
elseif size(ImIn,3)==3
    redChan=ImIn(:,:,1); redChan(BWoutline)=linecolor(1);
    greenChan=ImIn(:,:,2); greenChan(BWoutline)=linecolor(2);
    blueChan=ImIn(:,:,3); blueChan(BWoutline)=linecolor(3);
    outlined_mask(:,:,1)=redChan;
    outlined_mask(:,:,2)=greenChan;
    outlined_mask(:,:,3)=blueChan;
end

