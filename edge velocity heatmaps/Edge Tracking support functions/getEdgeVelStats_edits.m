function [protrusions,retractions,image]=getEdgeVelStats_edits(celldir, uniform_end,protThresh,retrThresh,sizeThresh)
% Input: 
%    start - end of uniform growth, after which symmetry breaking will occur 
%    uniform_end - when cells have stopped spreading
%    retrThres - edge velocity threshold for protrusions
%    sizeThres - size threshold for protrusions
% Output:
%     Image
%     protrusions, retractions: Structure arrays with the following parameers: 
%     - Area
%     - Location (centroid)
%     - Image of protrusion (logical, with smallest possible outline)
%     - PixelValues
%     - x - duration
%     - y - size in number of windows
%     - median velocity. 
%     - Image with protrusions/retractions outlined. 

% Arnold Hayer, July 26, 2017
%edited from original Seph Marshall, June 5, 2019

% retrThresh=-5;
% protThresh=5;

% cells = [1,2,4,5,8,9,10,11]; 


%% Filter out small junk
root= celldir;

load('CMAP_blue_grey_yellow.mat');
load([root,filesep,'Protrusion and FRET Values.mat'],'protvalsWindowF','fretvalsF'); 

edgemap=protvalsWindowF;
edgemapNew=protvalsWindowF(:,uniform_end:end);

fretmap=fretvalsF;
fretmapNew=fretvalsF(:,uniform_end:end);

protMap=edgemapNew>protThresh;
protMap=protMap.*edgemapNew;

protMapFull=edgemap>protThresh;
protMapFull=protMapFull.*edgemap;

retrMap=edgemapNew <retrThresh;
retrMap=retrMap.*edgemapNew;

retrMapFull=edgemap<retrThresh;
retrMapFull=retrMapFull.*edgemap;



protMap1=bwareaopen(protMap,sizeThresh);protLabel=bwlabel(protMap);
retrMap1=bwareaopen(retrMap,sizeThresh);retrLabel=bwlabel(retrMap);
% 
% f1 = figure; 
% 
% %ax1=subplot(1,2,1);
% imagesc(protMap,[-12 12]);
% title('Protrusions');
% colormap(f1,cmap); 
% 
% hold off; 
% f2 = figure; 
% % ax2=subplot(1,2,2);
% imagesc(retrMap,[-12,12]);
% title('Retractions');
% colormap(f2,cmap);
% hold off; 

%f3= figure;
%  subplot(1,2,1);imagesc(LOprot);
%  subplot(1,2,2);imagesc(LOretr);
%im1=ratio2RGB(edgemapNew,[-12 12]);
%colourmap(f3, cmp); 
%image=DrawMaskOutline(im1,retrMap1,[0 0 0]);
 %imagesc(image);
%colormap(image, cmap); 
% protrusions
protrusions=regionprops(protMap1,edgemapNew,'Area','Centroid','Image','PixelValues','PixelList','PixelIdxList','FilledImage');

% retractions
retractions=regionprops(retrMap1,edgemapNew,'Area','Centroid','Image','PixelValues','PixelList','PixelIdxList','FilledImage');

for k=1:size(protrusions,1)
   % imagesc(protrusions(k).Image);
    protrusions(k).coorRange=max(sum(protrusions(k).Image,1));
    protrusions(k).frameRange=max(sum(protrusions(k).Image,2));
    protrusions(k).mspeed=median(protrusions(k).PixelValues);
    protrusions(k).framestart=min(protrusions(k).PixelList(:,1))+uniform_end;
    protrusions(k).frameend=max(protrusions(k).PixelList(:,1))+ uniform_end;
    protrusions(k).coorstart=min(protrusions(k).PixelList(:,2));
    protrusions(k).coorend=max(protrusions(k).PixelList(:,2));
    protrusions(k).startref=uniform_end; 
    protrusions(k).numframes =length(protvalsWindowF); 
    

    FretTable=zeros(size(protrusions(k).PixelList,1),1);
  % MyosinTable=zeros(size(protrusions(k).PixelList,1),1);
   

    for i=1:size(protrusions(k).PixelList,1)
        intx=protrusions(k).PixelList(i,1);
        inty=protrusions(k).PixelList(i,2);
        
        FretTable(i,1)=fretvalsF(inty,intx+uniform_end);
     %  MyosinTable(i,1)=myosin(inty, intx+uniform_end);
         
    end
    
    s=size(protrusions(k).PixelList,1);
   
 
 
    protrusions(k).fretavg=mean(FretTable);
 %  protrusions(k).myosinavg=mean(MyosinTable);
  
end

% filter out events less than 5 frames, they tend to be mistakes anyways 
protrusions= protrusions(arrayfun(@(x) x.frameRange >= 5, protrusions));

save([root,filesep,'protrusionlist.mat'],'protrusions');

for k=1:size(retractions,1)
   % imagesc(retractions(k).Image);
    retractions(k).coorRange=max(sum(retractions(k).Image,1));
    retractions(k).frameRange=max(sum(retractions(k).Image,2));
    retractions(k).mspeed=median(retractions(k).PixelValues);
    retractions(k).framestart=min(retractions(k).PixelList(:,1))+uniform_end;
    retractions(k).frameend=max(retractions(k).PixelList(:,1))+ uniform_end;
    retractions(k).coorstart=min(retractions(k).PixelList(:,2));
    retractions(k).coorend=max(retractions(k).PixelList(:,2));
    retractions(k).startref=uniform_end; 
    retractions(k).numframes =length(protvalsWindowF); 
    
    
 FretTable2=zeros(size(retractions(k).PixelList,1),1);
  %MyosinTable2=zeros(size(retractions(k).PixelList,1),1);
 
    for i=1:size(retractions(k).PixelList,1)
        intx=retractions(k).PixelList(i,1);
        inty=retractions(k).PixelList(i,2);
        FretTable2(i,1)=fretvalsF(inty,intx+uniform_end);
  %     MyosinTable2(i,1)=myosin(inty,intx+uniform_end);
        
    end 
        
  
       s=size(retractions(k).PixelList,1);
   
 
    
    retractions(k).fretavg=mean(FretTable2);
  %  retractions(k).myosinavg=mean(MyosinTable2);
  
end

retractions= retractions(arrayfun(@(x) x.frameRange >= 5, retractions));
save([root,filesep,'retractionlist.mat'],'retractions');




%% mapping protrusions/retractions to location on edge velocity map 
% optional for visualization
% hold on; 
% axis tight; 
% axis ij;
% 
% imagesc(protMapFull,[-12,12]);
% title('Labelled Protrusions');
%  colormap(cmap);
%  
%   xline(uniform_end, 'w','END UNIFORM')
%     
% for j = 1:size(protrusions,1)
%    x=protrusions(j).Centroid(1,1);
%    y=protrusions(j).Centroid(1,2);
%     text(x+uniform_end, y, num2str(j),'color','b');
% end
% 
% hold off; 
% 
% hold on; 
% axis tight; 
% axis ij;
% 
% imagesc(retrMapFull,[-12,12]);
% title('Labelled Retractions'); 
%  colormap(cmap);
%  xline(uniform_end, 'w','END UNIFORM');
%     
% for w = 1:size(retractions,1)
%     
%    x=retractions(w).Centroid(1,1);
%    y=retractions(w).Centroid(1,2);
%     text(x+uniform_end, y, num2str(w),'color','w');
% end 
% 
% hold off; 
end 
% 
