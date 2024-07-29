clear; clc;

% Specify folder where raw background images are located:
bgpath='H:\180620_IXM\background'; 

% Specify the channels that were used, and which are part of the filenames
channels={'TRITC' 'GFP'};
%channels={'CFP' 'DAPI' 'mCherry'};

for chan=1:numel(channels);
    calculate_bg_img_rm_blobs(bgpath,channels{chan});
    disp(channels(chan));
    disp(num2str(chan));
end

disp('done!');