function [im_ch1, im_ch2] = nd2read_hayer_2chan(filename,finfo,counter,varargin)
tic
%finfo = nd2finfo(filename);
%disp(['analyzing file structure used ', sprintf('%0.2f', toc), ' seconds'])

c = numel(num2str(counter)); 

im_ch1 = zeros(finfo.img_width, finfo.img_height, 'uint16');
im_ch2 = zeros(finfo.img_width, finfo.img_height, 'uint16');


fid = fopen(filename, 'r');
fseek(fid, finfo.file_structure(strncmp(strcat('ImageDataSeq|',num2str(counter),'!'), ...
  {finfo.file_structure(:).nameAttribute}, 14+c)).dataStartPos, 'bof');

tic
% Image extracted from ND2 has image width defined by its first dimension.
if finfo.padding_style == 1
  
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch1(:, ii) = temp(1, :);
        im_ch2(:, ii) = temp(2, :);
    
        fseek(fid, 2, 'cof');
    end
end

  
    for ii = 1: finfo.img_height
        temp = reshape(fread(fid, finfo.ch_count * finfo.img_width, '*uint16'),...
          [finfo.ch_count finfo.img_width]);
        im_ch1(:, ii) = temp(1, :);
        im_ch2(:, ii) = temp(2, :);
   
    end 
 


fclose(fid);

im_ch1 = permute(im_ch1, [2 1]);
im_ch2 = permute(im_ch2, [2 1]);

  

%disp(['reading complete image data used ', sprintf('%0.2f', toc), ' seconds'])
end