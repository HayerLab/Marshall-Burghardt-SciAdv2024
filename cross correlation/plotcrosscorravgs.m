%Cross correlation Averages 

%% use this section to get the averages out of an entire set 
clc; clear; 
root='f:\example dataset ';
cells = [1:1:20]; %number of cropped cells in dataset 

cell_arr=cell(1,size(cells,2),1);

%%

depth =6; %specify edge depth data you wish to analyze 

for loop=1:size(cells,2)
   
    fileKey=strcat(num2str(cells(1,loop)));
    
    start =  2; 
    
    load([root,filesep,fileKey,filesep,'output', filesep, 'edge_vels', filesep, 'edge vel mapping_',num2str(depth),filesep,'Protrusion and FRET values.mat'],'ezrinF','protvalsWindowF','membranecytoF'); %'myosin', 'cytoF')
    
  
 %velocity vector created by averaging vectors before and after frame of
 %interest
 
    edgeVel_arr =protvalsWindowF;  
    temporary = NaN(180,size(protvalsWindowF,2)-1);
    for row = 1: size(protvalsWindowF,1)
        for col = 1:size(protvalsWindowF,2)-1
          temporary(row,col)=(protvalsWindowF(row,col)+protvalsWindowF(row,col+1))/2;
        end 
    end 
     edgeVel_arr_adjusted = temporary(:,start-1:end); 
     protExp_arr=efretvalsF(:,start:end-1);
    %because edge Velocity map now 2 frames smaller than protein expression (first and last frame cut off);
    
   protExp_arr(isnan(protExp_arr))=0; 
    
    
    %using CLT to shift each distribution to a normal distribution
    
    Z_edge=( (edgeVel_arr_adjusted-nanmean(edgeVel_arr_adjusted(:))) / std(edgeVel_arr_adjusted(:))); 
    Z_prot=( (protExp_arr-nanmean(protExp_arr(:))) / std(protExp_arr(:)));
     
    
    %generating cross correlation for each specific coordinate window 
    avg_table = NaN(180,size(protExp_arr,2)); 
     for i = 1:size(Z_edge,1) % 
  
     [c,lags] = xcorr((Z_edge(i,:)),Z_prot(i,:), 'coeff');
   
     avg_table(i,1:size(lags,2))=c; 
    
     end 
     
  %compiling average for a given cell  
  for c = 1:size(avg_table,2)
      meanval(1,c) = nanmean(avg_table(:,c)); 
  end    
      %discarding xcorr data outside the range of [-20,20] 
      for i= 1:size(lags,2) 
     
          if (abs(lags(1,i)) > 20)
              meanval(1,i)= NaN;  
   
          end 
          
      end 
      
      %now all mean vals should be centred around zero with range [-20,20] 
      meanval=meanval(~isnan(meanval));
     
      %inputting this into cell array that stores 
      cell_arr{1,loop}= meanval;
  
end 


   
%% plot all xcorrs from a trial 

overall_avg = zeros(1,41); 
for c =1:41
  
    counter =0;
    for cell =1:size(cell_arr,2)
        
 temp = cell_arr{1,cell};
        
  counter = counter + temp(1,c);   
    end 
    
    overall_avg(1,c) = counter/size(cell_arr,2); 
end 

f1=figure; 


hold on; 
%grid on; 
    title('Edge Vel vs.Rho FRET XCorr');
    xlabel('Lag (min)','FontWeight','bold');
    ylabel('Correlation Coefficient', 'FontWeight','bold'); 
    xline(0, '--'); 
    yline(0,'--');
   
for x= 1:size(cell_arr,2)
    
    plot([-20:20],cell_arr{1,x});
    
end 

plot([-20:20],overall_avg,'Color',[0,0,0], 'LineWidth', 3 ); 

yline(0,'--');
xline(0,'--'); 
xticks([-19.2 -16.8 -14.4 -12.0 -9.6 -7.2 -4.8 -2.4 0 2.4 4.8 7.2 9.6 12.0 14.4 16.8 19.2])
xticklabels({'-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8'}); 

hold off; 

% save data where desired 