clear;clc;
datadir='F:\170619_3i\data_170726';
load([datadir,filesep,'EdgedynamicsData_sheet_170619.mat']);
%%
%edgeData=edgeData(1:6); % for retraction stats
for cellID=1:size(edgeData,2)
    
%     currIm=ratio2RGB(edgeData(cellID).protvalsWindowF,[-12 12]);
%     imshow(currIm);
    currCellprotStats=edgeData(cellID).protrStats;
    currCellretrStats=edgeData(cellID).retrStats;
    
%     % Plot cumulative protrusion for each cell
%     %subplot(5,7,cellID);
%     protvalsWindowFcurr=cell2mat({edgeData(cellID).protvalsWindowF}');
%     protsum=sum(protvalsWindowFcurr,2);
%     polarplot((1:2:360)*pi/180,protsum);rlim([-300 300]); hold on;
    
     % for retractions
    meanProtDuration(cellID,1)=mean(cell2mat({currCellprotStats.xmax}'))/4; % frequency 4 frames/min
    meanProtSize(cellID,1)=mean(cell2mat({currCellprotStats.ymax}'))/180; % in % of cell circumference, 180 windows
    numProtr(cellID,1)=length(currCellprotStats)*3600/(90*15); % per h
    for protID=1:length(currCellprotStats)
        ProtrSpeed(protID,1)=mean(currCellprotStats(protID).PixelValues);
    end
    meanProtrSpeed(cellID,1)=mean(ProtrSpeed,1);
    
    meanRetrDuration(cellID,1)=mean(cell2mat({currCellretrStats.xmax}'))/4;
    meanRetrSize(cellID,1)=mean(cell2mat({currCellretrStats.ymax}'))/180;
    numRetr(cellID,1)=length(currCellretrStats)*3600/(90*15);
    for protID=1:length(currCellretrStats)
        RetrSpeed(protID,1)=mean(currCellretrStats(protID).PixelValues);
    end
    meanRetrSpeed(cellID,1)=mean(RetrSpeed,1);
end
disp('done!');

%%

