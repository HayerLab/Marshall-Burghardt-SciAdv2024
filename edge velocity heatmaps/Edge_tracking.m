%% Test protrusion/FRET measurements
% Required functions: 
%       functions located folder "trackingcode"  
%       parametrizeCellEdge.m
%       computeProtrusionValues.m
%       getWindowLabelMap.m       
%       mapSmoothedEdgeCoorsToCellBoundary.m
%       matrixCircularlyPermuteRows.m
%       computeNormalVectorsFromParametrizedCellEdge.m
%       CMAP_blue_grey_yellow.mat
%       ndnanfilter.m
%       vect.m


%need to add "Statistics and Machine Learning Toolbox" in MATLABS add-on
%section. 

%saves following figures: labelled cell mask with coordinates, area change,
%and edge velocity heat maps 
clc; clear; 
cells = [1];   %number of cropped cells that you wich to analyze 

for place=1:size(cells,2)
cells = [1];   
     


for depth = 1:5 % number of edge depths you wish to analyze 


 %for 20x no binning, or 40x 2x2 bin etc. 0.325 um/px
depths = [6; 



 
root='F:\example_dataset_240522';
rawdir=([root,filesep,'cropped',filesep, strcat(num2str(cells(1,place))),filesep,'output']); %
datadir=([rawdir,filesep,'edge_vels', filesep,  strcat('edge vel mapping_',num2str(depths(1, depth)))]); 
if ~exist(datadir)
    mkdir(datadir)
end 


% Parameters for cell edge parametrization
nFretWindows=180;                                   % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam=5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam=nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing=10;                                   % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge

edgeDepthDist=depths(1, depth);   % Number of pixels deep for the windows for computing FRET values.


binning=1;            %only change if binning is changed while using same objective magnification!!




% Load previously determined sequence of masks


load([rawdir,filesep,'RatioData_raw.mat']);
load([rawdir,filesep,'ezrin_data', filesep,'CytoRatioData.mat']);


%% %% delete unwanted frames (optional!!) 
% ie if cell touches another cell at the end of the movie and you wish to
% exclude it from edge tracking analysis  

% flashFrame=141:150; 
% maskFinal(flashFrame)=[];
% imCAAXOutline(flashFrame) = []; 
% ezrin_ratio(flashFrame)=[];
% membrane_cyto_ratio(flashFrame)=[];
% cellCoors(flashFrame)=[]; 
% % im_mRuby(flashFrame) = []; 


%% determines cell trajectories. 
% Will create multiple if there are multiple masked cells per crop 
% in your loadeddata 

%first checks if any trajectories exist, and then skips through any empty
%cells at the beginning 

%begins edge tracking at the first non-empty cell in the array cellCoors


empty_count=0;
if  ~isempty(cellCoors)

    for i=1:size(cellCoors,2)
       temp = cellCoors(1,i);
       contents=temp{1,1};
           if isempty(contents)
             empty_count=empty_count+1;
             continue; 

           else   
            traj=ultTrackAnnSearch(cellCoors(1,empty_count+1:end),'pairrule','fwdbckmtch','maxdisp',100,'minlength',5,'verbose',false); 
         end 
    end 
else
    traj={};
end
fprintf('%i trajectories.\n',length(traj));






%% Manually select cell for analysis 
selectedCell=1; %input which trajectory
isConnect = false; %true if you had to connect broken trajectories

thisTraj=traj{selectedCell};
centroid_coordinates=zeros(2,size(imCAAXOutline,2));
cell_area=0; 
if ~(isConnect) 
    start =thisTraj(1,5);
else
    
    start=thisTraj(1,4);

    for z = start:start+size(thisTraj,1) -1
        thisTraj(z,4)=z;
    end 

end

index = 1;
for imnum=start:start+size(thisTraj,1) -1
   
        % Make a one cell mask
        objects=regionprops(maskFinal{empty_count+imnum},'PixelIdxList','PixelList','Centroid','BoundingBox', 'Area');

        cellCent=round(thisTraj(find(thisTraj(:,end)==imnum,1),1:2));
        centroid_coordinates(:,imnum) = cellCent'; 
       
       if size(objects,1) == 1
        cell_area =cell_area +objects.Area; 
      end 

      
       if size(objects,1) > 1
        
           cell_areas = [objects.Area]; 
           cell_area = cell_area+max(cell_areas); 
           
       end 
           

            cellNum=find(round(arrayfun(@(x) x.Centroid(1),objects))==cellCent(1) & round(arrayfun(@(x) x.Centroid(2),objects))==cellCent(2));
            thisMask_raw=false(size(maskFinal{imnum}));
            thisMask_raw(objects(cellNum).PixelIdxList)=true;        
            thisMask_raw=false(size(maskFinal{imnum})); % sets a frame of the size of the mask to zero
            thisMask_raw(objects(cellNum).PixelIdxList)=true; % sets values of object 2 to 1

            thisMask = thisMask_raw;
            

                % Parametrize cell edge and compute protrusion values
                if index ==1
                    [edgeCoors{index}, edgeCoorsSmoothed{index}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning));
                else
                    [edgeCoors{index},edgeCoorsSmoothed{index}]=parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning),edgeCoors{index-1});
                   
                    protvals(:,index-1)=vect(computeProtrusionValues(edgeCoors{index-1},edgeCoorsSmoothed{index-1},edgeCoors{index}));
                end
                windowCoors{index}=edgeCoors{index}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);

            % Compute edge regions and local fret ratio values
            labelMask{index}=getWindowLabelMap(thisMask,windowCoors{index},(edgeDepthDist/binning));
         
            for k=1:size(windowCoors{index},1)
                
                %comment out as needed for specific data 
                 fretvals(k,index)=mean(imRatio_raw{index+empty_count}(labelMask{index}==k));
                 myosin(k,index)=mean(im_mRuby{index+empty_count}(labelMask{index}==k));
                 ezrin(k,index)=mean(ezrin_ratio{index+empty_count}(labelMask{index}==k));
                 membrane_cyto(k,index)=mean(membrane_cyto_ratio{index+empty_count}(labelMask{index}==k));
          
  
            end
           MaskThisCell{index}=thisMask;  
   
           disp(num2str(imnum));

           index=index+1;             
end

protvalsWindow=zeros(nFretWindows,size(protvals,2));

 for k=1:edgeOversamplingParam
     
     protvalsWindow=protvalsWindow+protvals(k:edgeOversamplingParam:end,:);
     
 end
 %% calculate total distance travelled by cell centroid and average cell area 
 distance = 0; 
 for i = 2:size(centroid_coordinates,2)
     dx = centroid_coordinates(1,i)-centroid_coordinates(1,i-1); 
     dy = centroid_coordinates(2,i)-centroid_coordinates(2,i-1); 
     
  distance = distance+ sqrt((dx^2)+(dy^2));       
 end 
 
 %% filtered protusionvalues
 
   % filtered protrusion values
protvalsWindowF=ndnanfilter(protvalsWindow,fspecial('disk',2),'replicate');
fretvalsF=ndnanfilter(fretvals,fspecial('disk',2),'replicate'); 
myosinF=ndnanfilter(myosin,fspecial('disk',2),'replicate');
ezrinF=ndnanfilter(ezrin,fspecial('disk',1),'replicate');
 membranecytoF=ndnanfilter(membrane_cyto,fspecial('disk',1),'replicate');

 
%% Plot maps - with thresholds
load('desired colormap.mat'); %any existing colormap can be used for visualization
load('CMAP_blue_grey_yellow.mat'); %default one used in paper 
close all;


f1=figure; 


hold on;


%following code creates 2 by 2 subplot where edge velocity, FRET activity, 
% and other protein localization can be observed. Must specify which data
% you want to plot where 

ax1=subplot(2,2,1);imagesc(protvalsWindow,[-13,13]);title('Edge Velocity');
colormap(ax1,BCWOR); 
%25 s intervals 
   xticks([0 24 48 71 95 120 144])
   xticklabels({'0','10','20','30','40','50' '60'});

ax2=subplot(2,2,2);imagesc(ezrinF,[-1 3] );title('Ezrin');
xlim([0 150]); 
 xticks([24 48 71 95 120 144])
  xticklabels({'0','10','20','30','40','50'});


ax3=subplot(2,2,3);imagesc(protvalsWindowF, [-13 13]); title ('Filtered');
colormap(ax3,BCWOR);
xlim([0 150]); 
xticks([ 24 48 71 95 120 144])
xticklabels({'0','10','20','30','40','50'});

 ax4=subplot(2,2,4);imagesc(myosinF,[0 3]); title ('Myosin localization');
 xlim([0 150]); 
 xticks([0 24 48 71 95 120 144])
 xticklabels({'0','10','20','30','40','50' '60'});

 hold off; 
saveas(f1,strcat(datadir,'\','edge_velocity_map.png'))
saveas(f1,strcat(datadir,'\','edge_velocity_mapM.fig'))


% save all new data into mat file 

save(strcat(datadir,'\','Protrusion and FRET Values.mat'),'protvalsWindow','protvalsWindowF','distance', 'ezrin', 'ezrinF', 'membrane_cyto','membranecytoF', 'fretvals', 'fretvalsF'); % 'avg_cell_area'%'myosin','myosinF'); % ;


%close all; clc;

end 
end 
 clear; clc; 

 %% count number of protrusions and retraction event/cell 
 % this section calls on function getEdgeVelStats_edits 
 % which identifies number of protrusion and retraction events per cell
 % that are above a certain size and speed threshold 

 % getEdgeVelStats_edits(celldir, startframe, prot_thresh, retr_thresh,
 % size thresh); 

 %celldir = filepath of data 
 %startframe: frame at which you wish to begin counting. Default one. 
 %prot_thresh: min avg speed in px/frame for protrusions 
 %retr_thresh: max avg speed in px/frame for retractions 
 %size thresh: minimum event size in pixels squared of event from
 %spatiotemporal heat map for event to be included S

 cells = [1:1:20]; 

 depths = [3]; %specify edge depth you wish to use 
 
 root = 'F:example_dataset_240522\cropped';
 for i = 1:20
    
     for j= 1:size(depths,2)
         
         
      celldir = ([root, filesep, strcat(num2str(i)),filesep, strcat('output', filesep,'edge_vels', filesep, 'edge vel mapping_', num2str(depths(1,1)))]);       
         
         
     getEdgeVelStats_edits(celldir,1,2.5,-2.5,25);
     
     
     end
     
 end 
 
