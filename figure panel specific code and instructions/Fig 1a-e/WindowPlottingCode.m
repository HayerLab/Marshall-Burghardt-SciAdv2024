% code used to plot edge coordinate windows over a cell image 


% Need to load a RatioData.mat file generated using singlecellprocessing.m

% Note: Make sure you have the 'Statistical Learning
               % Toolbox' by Dahua Lin installed. If not, you can add it by
               % HOME>Add-Ons>Manage Add-Ons>Get Add-Ons. 
               % Reference: Dahua Lin (2020). Statistical Learning Toolbox 
               % (https://www.mathworks.com/matlabcentral/fileexchange/12333-statistical-learning-toolbox), 
               % MATLAB Central File Exchange. Retrieved November 17, 2020.

stack = readTIFFstack("exampledatalocation\Rho-FRET.tif");

mask = zeros(size(maskFinal{1,1}));
for frame = 1:size(maskFinal,2)
    mask(:,:,frame) = maskFinal{1, frame};
end

%% Detect trajectories

%first checks if any trajectories exist, and then skips through any empty
%cells at the beginning 

%added iterations that deal with empty frames before cells are dropped,
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
           
               
            traj=ultTrackAnnSearch(cellCoors(1, empty_count+1:end),'pairrule','fwdbckmtch','maxdisp',100,'minlength',5,'verbose',false);
   
         end 
    end 
else
    traj={};
end
fprintf('%i trajectories.\n', length(traj));

%%

% Parameters for cell edge parametrization
nFretWindows = 180;                                   % Number of windows to use for FRET vs edge correlation measurements
edgeOversamplingParam = 5;                            % How many times more points should the edge of the cell be tracked at ("subwindow")
nPointsParam = nFretWindows*edgeOversamplingParam;    % Number of points to track on the cell edge
pdSmoothing = 10;                                     % Used with imclose to make the selection of points for tracking less dependent on noise or wrinkles in the cell edge


edgeDepthWidth = 29; % The distance in pixels between the outer and inner limits of the window band
edgeDepthDist = 0;% The distance in pixels between the cell edge and the outer limit of the window band

startFrame=1;
endFrame=130;
binning = 1;            %only change if binning is changed while using same objective magnification!!


%% Manually select cell for analysis 

selectedCell = 1; %input which trajectory
isConnect = true; % if you had to connect broken trajectories
% deleteFrame = 46; % manually input frame that needs to be deleted
thisTraj = traj{selectedCell};

if ~isConnect %modifying this Traj code if you needed to correct broken trajectores
    start = thisTraj(1,6);
else
    
    start = thisTraj(1,5);

    for z = start:start+size(thisTraj,1) -1
        thisTraj(z,5)=z;
    end 
end

%%
index = 1;
for imnum=start:start+size(thisTraj,1) -1
   
    
        objects=regionprops(maskFinal{empty_count+imnum},'PixelIdxList','PixelList','Centroid','BoundingBox');    

        cellCent=round(thisTraj(find(thisTraj(:,end)==imnum,1),1:2));


            cellNum = 1;
            thisMask_raw=false(size(maskFinal{imnum})); % sets a frame of the size of the mask to zero
            thisMask_raw(objects(cellNum).PixelIdxList)=true; % sets values of object 2 to 1

            % Smooth mask
            thisMask=imerode(thisMask_raw,strel('disk',4));
            thisMask=imdilate(thisMask,strel('disk',4));
            thisMask=bwareaopen(thisMask,300);

                % Parametrize cell edge and compute protrusion values
                if index == 1
                    [edgeCoors{index}, edgeCoorsSmoothed{index}] = parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning));
                
                else
                    [edgeCoors{index}, edgeCoorsSmoothed{index}] = parametrizeCellEdge(thisMask,nPointsParam,round(pdSmoothing/binning),edgeCoors{index-1});
                    %windowCoors{frameCount}=edgeCoors{frameCount}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);
                    [protrusionValues, protrusionVectors] = computeProtrusionValues(edgeCoors{index-1}, edgeCoorsSmoothed{index-1}, edgeCoors{index}, edgeCoorsSmoothed{index});  
                    protvals(:,index-1)=vect(protrusionValues);
                    prVect{index} = protrusionVectors;
%                     normV{index} = computeNormalVectorsFromParametrizedCellEdge(edgeCoorsSmoothed{index-1});
                end
                normV{index} = computeNormalVectorsFromParametrizedCellEdge(edgeCoorsSmoothed{index});
                windowCoors{index}=edgeCoors{index}((edgeOversamplingParam+1)/2:edgeOversamplingParam:end,:);

            % Compute edge regions and local fret ratio values
            labelMask{index}=getWindowLabelMap(thisMask,windowCoors{index},(edgeDepthWidth/binning),edgeDepthDist); %
         
%             for i=1:size(windowCoors{index},1)
                
%                 fretvals(i,index)=mean(imRatio{index+empty_count}(labelMask{index}==i));
%                 myosin(k,index)=mean(im_mRuby{index+empty_count}(labelMask{index}==k));
%                 actin(k,index) = mean(im_actin{index + empty_count}(labelMask{index} == k));
%                 actin(k,index) = mean(masked_stack{index + empty_count}(labelMask{index} == k));
          
%             end
           MaskThisCell{index}=thisMask;  
  

           index=index+1;             
end

protvalsWindow=zeros(nFretWindows,size(protvals,2));

 for i=1:edgeOversamplingParam
     
     protvalsWindow=protvalsWindow+protvals(i:edgeOversamplingParam:end,:);
     
 end
 
 %% Select a frame
 
figure; 
imagesc(protvals, [-15,15]); colorbar;
caxis([-2,2]); axis square;
            

for frameNum=1
     close all; 
    %imagesc(MaskThisCell{frameNum}); hold on;
    
    image = imFRETOutline{1,frameNum};
   %image = immRuby_outline{1,frameNum};
   % image =imRatio{1, 1};   
   %image=ratio2RGB( image,[0.7 1.3]);
      h = figure('visible','on');

      axis ij; 
    imagesc(image);
    imshow(image);
 
    
 
    hold on; 
    
     rectangle('Position', [380 55 216 216],'EdgeColor', [1 1 1]); 
     rectangle('Position', [390 65 45 5],'EdgeColor', [1 1 1], 'FaceColor', [1 1 1]); 
%     
 for frame = frameNum
    for w = 1:nFretWindows
        windowMask = (labelMask{1, frame} == w);        
        [r1,c1] = ind2sub(size(windowMask), find(windowMask, 1));
        temp = bwtraceboundary(windowMask, [r1,c1], 'N');
        plot(temp(:,2), temp(:,1), 'r', 'LineWidth', 0.5);
        windowOutline{frame, w} = temp;
    end
end
     
    normV=computeNormalVectorsFromParametrizedCellEdge(edgeCoorsSmoothed{frameNum});
    protvals=(computeProtrusionValues(edgeCoors{frameNum},edgeCoorsSmoothed{frameNum},edgeCoors{frameNum+1}, edgeCoorsSmoothed{frameNum+1}));
    %plot(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1)); quiver(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1),0.2*normV(:,2),0.2*normV(:,1));
    protvect=normV.*protvals;
    
    % filter protvect with moving average of windowSize 5
    windowSize=5;
    b=(1/windowSize)*ones(1,windowSize);
    a=1;
    
    protvectF(:,1)=filter(b,a,protvect(:,1));
    protvectF(:,2)=filter(b,a,protvect(:,2));
    
    protvectF_binned=zeros(180,2); 
    for i = 1:180
        protvectF_binned(i,1) = mean(protvectF((i-1)*5+1:5*i, 1));   
        protvectF_binned(i,2)= mean(protvectF((i-1)*5+1:5*i,2)); 
    end 
    
     plot(windowCoors{frameNum}(:,2),windowCoors{frameNum}(:,1)); quiver(windowCoors{frameNum}(:,2),windowCoors{frameNum}(:,1),0.4*protvectF_binned(:,2),0.4*protvectF_binned(:,1),'Color',[0 1 0]);
    
    % adding this to try to make the vectors different colours 
    for k = 1:180
        if protvalsWindow(k, frameNum) >=5
          plot(windowCoors{frameNum}(k,2),windowCoors{frameNum}(k,1)); quiver(windowCoors{frameNum}(k,2),windowCoors{frameNum}(k,1),protvectF_binned(k,2),protvectF_binned(k,1),'Color',[1 1 0],  'LineWidth', 2, 'AutoScale','off');
       
        elseif protvalsWindow(k, frameNum)<= -5
          plot(windowCoors{frameNum}(k,2),windowCoors{frameNum}(k,1)); quiver(windowCoors{frameNum}(k,2),windowCoors{frameNum}(k,1),protvectF_binned(k,2),protvectF_binned(k,1),'Color',[0 0 1], 'LineWidth', 2,'AutoScale','off');
      
        else
   
         plot(windowCoors{frameNum}(k,2),windowCoors{frameNum}(k,1)); quiver(windowCoors{frameNum}(k,2),windowCoors{frameNum}(k,1),protvectF_binned(k,2),protvectF_binned(k,1),'Color',[1 1 1],'LineWidth', 2, 'AutoScale','off');
                
       end 
    end
    
    %plot binned vectors 
 
 
    %plot all vectors, not smoothened
     %  plot(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1)); quiver(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1),0.2*protvect(:,2),0.2*protvect(:,1),'Color',[1 0 0]);
 
  %plot all vectors, smoothened
    plot(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1)); quiver(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1),0.4*protvectF(:,2),0.4*protvectF(:,1),'Color',[0 1 0]);
  
   
    
    
    plot(windowCoors{frameNum}(:,2),windowCoors{1}(:,1));hold on;
    plot(edgeCoors{frameNum}(:,2),edgeCoors{frameNum}(:,1));
    plot(edgeCoorsSmoothed{frameNum}(:,2),edgeCoorsSmoothed{frameNum}(:,1));
    qp.edgePos=windowCoors{frameNum};
    qp.edgeVel=6*(windowCoors{frameNum+1}-windowCoors{frameNum});
    quiver(qp.edgePos(:,2)',qp.edgePos(:,1)',qp.edgeVel(:,2)',qp.edgeVel(:,1)');hold off;%,'Color',[1 1 1]
    pause;  

rectangle('Position', [185 124 25 5],'EdgeColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', 1); 
rectangle('Position', [185 124 5 5],'EdgeColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth', 1); 
rectangle('Position', [185 124 10 5],'EdgeColor', [1 1 1], 'EdgeColor', [0 0 0], 'LineWidth',1); 
rectangle('Position', [185 124 15 5],'EdgeColor', [1 1 1],'EdgeColor', [0 0 0], 'LineWidth', 1); 
rectangle('Position', [185 124 19 5],'EdgeColor', [1 1 1],'EdgeColor', [0 0 0], 'LineWidth', 1); 
rectangle('Position', [185 124 22 5],'EdgeColor', [1 1 1],'EdgeColor', [0 0 0], 'LineWidth', 1);

 set(h, 'Renderer', 'painters');
% ax = gca; 
% 
% ax.YLim= [18 143]; 
% ax.XLim= [135 235];

cropped_graph=imcrop(image,[135, 18, 100, 125]) 
print(gcf,'-vector','-dsvg',['C:\Users\marsh\OneDrive - McGill University\research paper\results good_Feb2023\Fig 1\fullcell_wvectors_take2.svg']); % svg
end


                                       
%% Plot window traces
for frame = 80
    for w = 1:nFretWindows
        windowMask = (labelMask{1, frame} == w);        
        [r1,c1] = ind2sub(size(windowMask), find(windowMask, 1));
        temp = bwtraceboundary(windowMask, [r1,c1], 'N');
        plot(temp(:,2), temp(:,1), 'r', 'LineWidth', 2);
        windowOutline{frame, w} = temp;
    end
end