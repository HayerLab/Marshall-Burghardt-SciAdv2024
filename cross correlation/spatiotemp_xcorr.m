%% 2D cross correlation 
% creates single cell cross-correlation maps using the parula colourscale. 
% requires matlab files generated from the EdgeTracking_edits script. 

clear; 
clc; 

load('D:\example dataset\4\edge_vels\edge vel mapping_6\Protrusion and FRET Values.mat'); 

colour_array =jet(100); 
edgeVel_arr =protvalsWindowF;  
    temporary = NaN(180,size(protvalsWindowF,2)-1);
    for row = 1: size(protvalsWindowF,1)
        for col = 1:size(protvalsWindowF,2)-1
          temporary(row,col)=(protvalsWindowF(row,col)+protvalsWindowF(row,col+1))/2;
        end 
    end 
    edgeVel_arr_adjusted = temporary(:,1:end); 
    protExp_arr=fretvalsF(:,2:end-1); %because edge Velocity map now 2 frames smaller than protein expression (first and last frame cut off);
    
    protExp_arr(isnan(protExp_arr))=0; 
    
    
    %using CLT to shift each distribution to a normal distribution
    Z_edge=( (edgeVel_arr_adjusted-nanmean(edgeVel_arr_adjusted(:))) / std(edgeVel_arr_adjusted(:)));
   
     Z_prot=( (protExp_arr-nanmean(protExp_arr(:))) / std(protExp_arr(:)));
     
    avg_table = NaN(180,size(fretvalsF,2)); 
     for i = 1:size(Z_edge,1)
  
     [c,lags] = xcorr((Z_edge(i,:)), Z_prot(i,:),'coeff');
   
     avg_table(i,1:size(lags,2))=c; 
    
     end
     
     metric = size(lags,2)/2+0.5; 
 
  f= figure;  

    imagesc(avg_table(:, metric-20:metric+20),[-1,1]); 
    %colormap(f, colour_array); 
   % imagesc(f,[0 1]); 
    colormap(f); 
    xticks([0.8 3.2 5.6 8 10.4 12.8 15.2 17.6 20 22.4 24.8 27.2 29.6 32 34.4 36.8 39.2])
xticklabels({'-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8'}); 

xline(20,'--w','LineWidth',2);

xlabel('Lag (min)'); 
ylabel('Coordinate Window');

corr_table = avg_table(:, metric-20:metric+20); 
corr_avg = mean(corr_table); 

f1 = figure; 
 
plot([-20:1:20], corr_avg,'LineWidth',2);
xticks([-19.2 -16.8 -14.4 -12.0 -9.6 -7.2 -4.8 -2.4 0 2.4 4.8 7.2 9.6 12.0 14.4 16.8 19.2])
xticklabels({'-8' '' '-6' '' '-4' '' '-2' '' '0' '' '2' '' '4' '' '6' '' '8'}); 
xlabel('Lag (min)'); 
ylabel('Cross Correlation'); 

xline(0, '--k'); 
