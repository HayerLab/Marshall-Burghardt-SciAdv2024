function [ TiffStackOut] = readTIFFstack(TiffStackLocation)
% Speeds up reading of tif stacks
% copied from http://www.matlabtips.com/how-to-load-tiff-stacks-fast-really-fast/
% 
% Requires tifflib.mexmaci64 (for Mac) or tifflib.mexw64 (for Windows) to be within the search path
% 
% Input: Multiplane tiffstack, the type FIJI generates
% 
% Output: 3D array  
%
% 29 July 2015 Arnold Hayer

FileTif=TiffStackLocation;
InfoImage=imfinfo(FileTif);
mImage=InfoImage(1).Width;
nImage=InfoImage(1).Height;
NumberImages=length(InfoImage);
FinalImage=zeros(nImage,mImage,NumberImages,'uint16');
 
TifLink = Tiff(FileTif, 'r');
for i=1:NumberImages
   TifLink.setDirectory(i);
   TiffStackOut(:,:,i)=TifLink.read();
   %disp(num2str(i));
end
TifLink.close();
end

