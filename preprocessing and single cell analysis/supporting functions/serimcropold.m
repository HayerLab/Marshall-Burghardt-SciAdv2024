function [crop,rect] = serimcropold(series, ref)
% Crops a region of an image series, and outputs a smaller 3D matrix
% While loop ensures user is satisfied with choice
% If prime is 1, function will not accept a region with a prime dimension

if nargin < 2
    ref = series;
end
ok = 'n';
%figure;

while ok == 'n'
  %  colormap(hot(256));
   % imshow((double(mean(ref,3))))
    imshow(ref, [0 2000])
    axis image;
    title('Please select the region of interest','FontSize',18);
    [~, rect] = imcrop;
    close;
    rect = floor(rect);

    crop = uint16(zeros(rect(4)+1,rect(3)+1,size(series,3))); % Pre-allocates matrix

        for i=1:size(series,3)
            crop(:,:,i) = imcrop(series(:,:,i),rect);
        end


    ok = 'y';  
    % disp('You have selected the following area')
    % imgmovie(crop,0)
    % ok = input('Are you satisfied? (y/n)', 's');

    % if prime == 1
    %     if or(isprime(size(crop,1)),isprime(size(crop,2)))
    %     errordlg('You have selected an area with a prime dimension.  Please select again.','Fourier Transform Error')
    %     ok = 'n';
    %     end
    % end

end
