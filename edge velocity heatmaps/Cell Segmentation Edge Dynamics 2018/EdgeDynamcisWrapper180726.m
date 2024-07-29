%% Data input for 170619
% Cells were manually selected so as to only include cells without masking
% artifacts.

clear;
clc;


% Enter the well IDs ("shots") for each movie to be included in the analysis
shots={...
    '3_3_1'
        };

 % For each well ID ("shot"), add the cell IDs that should be included.    
 shotCells={
     [3 5]
        };


%% Parameters
addpath('E:\data_D\Stanford\Matlab\edge_tracking_RhoA2GCAAX\trackingcode');
dataIn='F:\180712_IXM';
dataOut='F:\180712_IXM';
if ~exist(dataOut)
    mkdir(dataOut);
end
   imYFPcrop
%% Parameters for cell edge parametrization
nFretWindows=180;                                   % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=10; %was 10                             % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge
edgeDepthDist=6; %was 10                            % Number of pixels deep for the windows for computing FRET values.

startFrame=1;
endFrame=90;
binning=1;
%% Cells in sheets
% Movies
cellCount=0;
for shotID=1:size(shots,1)
    %position=['CC-Rac1_sheet_15s_',num2str(single(shotID)),'_CFP-YFP-FRET'];
    load([dataIn,filesep,shots{shotID},'_MaskingData.mat']);
    % Tracking
    if min(cellfun(@(x) size(x,1),cellCoors))>0
    traj=ultTrackAnnSearch(cellCoors,'pairrule','fwdbckmtch','maxdisp',100,'verbose',false);
    else
        traj={};
    end
 % Cells   
    
    for cellID=1:size(shotCells{shotID},2)
        cellCount=cellCount+1;
        edgeData(cellCount).shotID=shots{shotID};
        edgeData(cellCount).cellID=cellID;
        thisTraj=traj{shotCells{shotID}(cellID)};
        % Make a one cell mask        
        frameCount=0;
        % Frames
        for imnum=startFrame:endFrame
           
            frameCount=frameCount+1;
            
            objects=regionprops(maskFinal{imnum},'PixelIdxList','PixelList','Centroid','BoundingBox');
            cellCent=round(thisTraj(find(thisTraj(:,end)==imnum,1),1:2));
            cellNum=find(round(arrayfun(@(x) x.Centroid(1),objects))==cellCent(1) & round(arrayfun(@(x) x.Centroid(2),objects))==cellCent(2));
            thisMask_raw=false(size(maskFinal{imnum}));
            thisMask_raw(objects(cellNum).PixelIdxList)=true;        
            thisMask_raw=false(size(maskFinal{imnum})); % sets a frame of the size of the mask to zero
            thisMask_raw(objects(cellNum).PixelIdxList)=true; % sets values of object 2 to 1

            % Smooth mask
            thisMask=imerode(thisMask_raw,strel('disk',4));
            thisMask=imdilate(thisMask,strel('disk',4));
            thisMask=bwareaopen(thisMask,500);% Clean up smoothing artifacts (remove objects smaller than 500 pixels). 

            % Parametrize cell edge and compute protrusion values
            if frameCount==1
                [edgeCoors{frameCount} edgeCoorsSmoothed{frameCount}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning));

                plot(edgeCoors{1,1}(:,1),edgeCoors{1,1}(:,2)); hold on; 
                plot(edgeCoorsSmoothed{1,1}(:,1),edgeCoorsSmoothed{1,1}(:,2))
            else
                [edgeCoors{frameCount},edgeCoorsSmoothed{frameCount}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning),edgeCoors{frameCount-1});
                %windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
                protvals(:,frameCount-1)=vect(computeProtrusionValues(edgeCoors{frameCount-1},edgeCoorsSmoothed{frameCount-1},edgeCoors{frameCount}));
            end
            windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);

            % Compute edge regions and local fret ratio values
            labelMask{frameCount}=getWindowLabelMap(thisMask,windowCoors{frameCount},edgeDepthDist/binning);
            YFPframe=imYFPcrop(:,:,frameCount);
            for k=1:size(windowCoors{frameCount},1)
                %fretvals(k,frameCount)=sum(imFRET(labelMask==k))/sum(imCFP(labelMask==k));
                fretvals(k,frameCount)=mean(imRatio_raw{frameCount}(labelMask{frameCount}==k));
                yfpvals(k,frameCount)=mean(YFPframe(labelMask{frameCount}==k));
            end
            MaskThisCell{frameCount}=thisMask;
            MaskThisCellRaw{frameCount}=thisMask_raw;
            disp(num2str(imnum));
        end % Frames

        protvalsWindow=zeros(size(fretvals)-[0 1]);
        for k=1:edgeOversamplingParam
            protvalsWindow=protvalsWindow+protvals(k:edgeOversamplingParam:end,:);
        end
         % Normalize protrusion data by dividing by the edgeOversampling factor and dividing by the time between images.
        protvalsWindow=(protvalsWindow/edgeOversamplingParam)./repmat(0.25*ones((endFrame-startFrame),1)',[nFretWindows 1]); % one frame every 15s, unit will be per minute 
        
        protvalsWindowF=ndnanfilter(protvalsWindow,fspecial('disk',2),'replicate');% filters the map
        % Normalize protrusion data by dividing by the edgeOversampling factor and dividing by the time between images.
        % protvalsWindow=(protvalsWindow/edgeOversamplingParam)./repmat(diff(timeData(startFrame:endFrame))',[nFretWindows 1]);
        
        edgeData(cellCount).MaskThisCellRaw=MaskThisCellRaw;
        edgeData(cellCount).MaskThisCell=MaskThisCell;
        edgeData(cellCount).protvalsWindow=protvalsWindow;
        edgeData(cellCount).protvalsWindowF=protvalsWindowF;
        edgeData(cellCount).fretvals=fretvals;
        edgeData(cellCount).yfpvals=yfpvals;
        
        % Analyze protrusions in each cell:
        [protrStats,retrStats,edgeDynImage]=getEdgeVelStats(protvalsWindowF,5,-3,10);
        edgeData(cellCount).protrStats=protrStats;
        edgeData(cellCount).retrStats=retrStats;
        imwrite(edgeDynImage,[dataOut,filesep,shots{shotID},'_Cell_#',num2str(shotCells{shotID}(cellID)),'_##',num2str(cellCount),'.tif'],'Compression','none');
        
    end % cells
    
  end % movies
save([dataOut,filesep,'EdgedynamicsData_180723.mat'],'edgeData','-v7.3');  

