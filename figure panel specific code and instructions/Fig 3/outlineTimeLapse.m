%root ='F:\Seph\research paper\Fig 1'; 

colour_array =jet(26);
 
 f1 = figure('visible','on');
      hold on; 
       axis ij;
 
 %change these depending on actual size of image - so it lines up perfectly with 
 % the number of pixels
 ylim([0 342]) 
xlim([0 456]);

%timepoints that you wish to plot (ie 1:1:150 for the whole video)
   for t = [51,58,61,63,65,68,75]
  
    temp = edgeCoors{1,t};
   plot(temp(:,2), temp(:,1),'Color',colour_array(t-50,:))
 
   end
   
   set(gca,'XColor', 'none','YColor','none')
hold off; 
  %  title('T = 28 min'); 
    
    
% CAN SEE ACTUAL CELL OUTLINE IF NEEDED HERE     
  f2=figure;  
  
 
  
         
 ylim([0 342])
xlim([0 456]);
  axis ij; 
image = imFRETOutline{1,51};
imagesc(image);
hold off; 

 
     