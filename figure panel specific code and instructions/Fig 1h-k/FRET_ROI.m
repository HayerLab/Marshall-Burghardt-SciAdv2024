%% FRET-ROI
% Interactive 

clear;clc;close all; warning off;

% Select on of the three following data folders.
root = '..\170801_3i';
%root='..\170810_3i';
%root='..\170920_3i';

rawdir=[root,filesep,'raw'];
datadir=[root,filesep,'data'];

files=getFilenames(rawdir,'_CFP Wide');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID=files{1}; % change index to loop through files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fileID=fileID(1:end-14);

CFPStack=double(readTIFFstack([rawdir,filesep,fileID,'_CFP Wide.tiff']));
FRETStack=double(readTIFFstack([rawdir,filesep,fileID,'_FRET Wide.tiff']));
REDStack=double(readTIFFstack([rawdir,filesep,fileID,'_Cherry Wide.tiff']));

CFPStack=ndnanfilter(CFPStack,fspecial('disk',2),'replicate');
FRETStack=ndnanfilter(FRETStack,fspecial('disk',2),'replicate');
REDStack=ndnanfilter(REDStack,fspecial('disk',2),'replicate');

RATIOStack=FRETStack./CFPStack;

REDgreyRange = [round(prctile(vect(REDStack(:,:,1)),3),1),round(prctile(vect(REDStack(:,:,1)),99),1)];
FRETgreyRange=[round(prctile(vect(FRETStack(:,:,1)),3),1),round(prctile(vect(FRETStack(:,:,1)),99),1)];

RGBim(:,:,1)=mat2gray(REDStack(:,:,1),REDgreyRange);
RGBim(:,:,2)=mat2gray(FRETStack(:,:,1),FRETgreyRange);
RGBim(:,:,3)=zeros(size(RGBim(:,:,1)));

%% Loop through control cells
i=[];
cellNum=1;
while isempty(i)
    if cellNum>1
        reply=input('NEXT CELL [ENTER], DONE [1] ');
        i=reply;
    end
    if isempty(i)
        subplot(1,2,1);
        imshow(mat2gray(REDStack(:,:,1),REDgreyRange));title(fileID);
        subplot(1,2,2);
        h_im=imshow(RGBim);title('Draw ROI on CONTROL cells, press [ENTER] to proceed to next cell or  [1/ENTER] to finish');    
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        ROI=imellipse;
        mask=createMask(ROI,h_im);
        for frameNum=1:size(RATIOStack,3)
            tempframe=RATIOStack(:,:,frameNum);
            ROIcontrol(frameNum,cellNum)=mean(vect(tempframe(mask)));
        end
        cellNum=cellNum+1;
    end
end
close all;

%% Loop through cells expressing both markers

i=[];
cellNum=1;
while isempty(i)
    if cellNum>1
        reply=input('NEXT CELL [ENTER], DONE [1] ');
        i=reply;
    end
    if isempty(i)
        subplot(1,2,1);
        imshow(mat2gray(REDStack(:,:,1),REDgreyRange));title(fileID);
        subplot(1,2,2);
        h_im=imshow(RGBim);title('Draw ROI on TEST cells, press [ENTER] to proceed to next cell or  [1/ENTER] to finish');    
        set(gcf, 'Position', get(0,'Screensize')); % Maximize figure.
        ROI=imellipse;
        mask=createMask(ROI,h_im);
        for frameNum=1:size(RATIOStack,3)
            tempframe=RATIOStack(:,:,frameNum);
            ROItest(frameNum,cellNum)=mean(vect(tempframe(mask)));
        end
        cellNum=cellNum+1;
    end
end
save([datadir,filesep,fileID,'_kinetics.mat'],'ROItest','ROIcontrol');

%% Plot results (only results of current timelapse are shown)
close all; 
subplot(1,2,1);plot(-135:15:600,ROIcontrol);hold on;
hold on;plot([0 0],[0 2]);
xlim([-150 600]);ylim([0.5 2]);title('control');

subplot(1,2,2);plot(-135:15:600,ROItest);
hold on;plot([0 0],[0 2]);
xlim([-150 600]);ylim([0.5 2]);title('test');