clear; clc; 
  %  uifig = uifigure(); 
     
    root = 'E:\seph backup\LOK SLK ERM KD';
    
    
    
    tiffDir = ([root, filesep,'tiff_stacks']); 
    
    bgCFP = ([ root, filesep, 'background\AVG_bgCFP.tif']); 
    bgFRET = ([root, filesep, 'background\AVG_bgFRET.tif']);
    
    cell= 1;
    filekey = '1_1_14'; 
    
     FRET = ([tiffDir,filesep, strcat(filekey,'_FRET_stacked.tif')]);  
     CFP= ([tiffDir,filesep, strcat(filekey,'_CFP_stacked.tif')]); 
     
     cellDir = ([root, filesep, 'cropped',filesep,'siERM', num2str(cell)]); 
     if ~exist(cellDir)
         mkdir(cellDir); 
     end 
    
       
            cropSite = 0;

            
            %FRET_stack = readFileToStack(FRET);
            
         FRET_stack=double(readTIFFstack(FRET));
         
         FRET_stack_1 = FRET_stack(:,:,1); 
          
         
         
         fg = figure;
         
%          imshow(FRET_stack_1,[0 2000]); 
%          
%          
         axis image;  %fg.WindowState = "maximized";
            

                
            
                [stackFRET, cropArea] = serimcropold(FRET_stack,mean(FRET_stack,3));
      
                
                Stack2TIFF(stackFRET, ['FRET.tif']);
                    
                    stackCFP = readFileToStack(CFP); 
                    stackCFP = imcrop3(stackCFP, [cropArea(1), cropArea(2), 1,...
                    cropArea(3), cropArea(4), size(stack,3)-1]);                    
                    Stack2TIFF(stackCFP, ['CFP.tif']);
            
                    
            close(fg);
            clear stack;
  
  
    close(uifig);
