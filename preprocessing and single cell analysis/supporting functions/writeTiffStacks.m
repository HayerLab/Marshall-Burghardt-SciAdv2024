function [] = writeTiffStacks(root, channels,row, col, numsites,timeFrames)
%writes raw organized microscopy data into TIFF stacks -- can be done
%manually using ImageJ also. 


%Seph Marshall 10/06/19

% inputs:
%root: string directory where raw tiff. images folder is found
%channels: string array, ie ["CFP" ; "FRET"];
%row: specific row from file tag 
%channel: specific channel from filetag 
%numsites: number of sites per well
%timeFrames: total number of frames in the experiment

%  outputs: 
% tiff stack to data directory specified below 

rawdir=[root,filesep,'raw'];
datadir=[root,filesep,'tiff_stacks'];

if ~exist(datadir)
    mkdir(datadir);
end


for s = 1:size(channels,1)
    
for i =1:numsites

     for j = 1:timeFrames

         dataID=strcat(num2str(row),'_',num2str(col),'_',num2str(i),'_',channels(s));
         dataKey=strcat(num2str(row),'_',num2str(col),'_',num2str(i),'_',channels(s),'_',num2str(j),'.tif');

      files=getFilenames(rawdir);
      files=files(boolRegExp(files, dataKey));
      files=convertStringsToChars(files);
        
        tempFrame=(imread([rawdir,filesep,files{1,1}])); 

       imwrite(tempFrame,[datadir,filesep,convertStringsToChars(dataID),'_stacked.tif'],'WriteMode','append','Compression','none');

     end   
    disp(num2str(i));

end

end 

end 
