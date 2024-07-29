%% Objects & Properties
for r = 1 % Fake for loop to enable section collapse

    ExpSettings.ObjectsPath = [ExpSettings.RootPath, 'Objects', filesep];
    if ~exist(ExpSettings.ObjectsPath, 'dir'); mkdir(ExpSettings.ObjectsPath); end    
    
    uifig = uifigure; uifig.WindowStyle = 'alwaysontop'; uifig.Name = 'Big cells, you are beautiful! ';
    
    for file = ExpSettings.Selected.files_to_analyze
        
        fileKey = getFileKey(file, 0, 0, metadata, ExpSettings.Selected.Convention);
        if logical(exist([ExpSettings.ObjectsPath, fileKey, '-Objects.mat'], 'file'))
            load([ExpSettings.ObjectsPath, fileKey, '-Objects.mat'], 'Objects');
            if length(Objects) ~= length(Options.Objects)
                disp('Error! Objects and Options.Objects should have the same size!');
                break
            end
        else % Initiate Objects with the right field names
            if exist('Objects', 'var'); clear Objects; end
            totalNumObjects = length(Options.Objects);
            Objects(totalNumObjects) = struct(); %#ok<SAGROW>
            for p = 1:totalNumObjects
                Objects(p).ExpID = ExpSettings.ExpID;
                Objects(p).File = 0;
                Objects(p).Site = 0;
                Objects(p).Crop = 0;
                Objects(p).FileKey = '';
                Objects(p).Treatment = '';
                Objects(p).Duration = '';
                Objects(p).DateMod = datestr(now);
                Objects(p).Frame = 1;
                Objects(p).ObjectIndex = 0;
                Objects(p).NumObjects2D = 0;
                Objects(p).CentroidX = NaN;
                Objects(p).CentroidY = NaN;
                Objects(p).Area = NaN;
                Objects(p).Perimeter = NaN;
                Objects(p).Solidity = NaN;
                Objects(p).PerSqrdOverArea = NaN;
                Objects(p).BBX = 1;
                Objects(p).BBY = 1;
                Objects(p).BBW = 1;
                Objects(p).BBH = 1;
                Objects(p).Eccentricity = NaN;
                Objects(p).Orientation = NaN;
                Objects(p).MajorAxis = NaN;
                Objects(p).MinorAxis = NaN;
                regProps = regionprops((1), 'PixelIdxList');
                Objects(p).PixelIdxList = {regProps.PixelIdxList};
                Objects(p).MeanIntensity = NaN;                  
            end            
        end
        
        if ~isempty(fieldnames(Objects)) && length(Objects) > 1
            % Sort Objects based on File, then Site, then Crop
            sTbl = struct2table(Objects);
            [tbl, ~]  = sortrows(sTbl, {'File', 'Site', 'Crop'});
            Objects = table2struct(tbl); clear index tbl sTbl;
        end


        index = 1;
        
        selectedSites = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file));
        indexes = [];
        for site = selectedSites
            indexes = [indexes, getObjectIndex(file, site, 0, Options.Objects)]; %#ok<AGROW>
        end
        
        if exist('Meta', 'var'); clear Meta; end
%         Meta(length(indexes)) = struct(); %#ok<SAGROW>
        if exist('tblObjects', 'var'); clear tblObjects; end
%         tblObjects(length(indexes)) = struct(); %#ok<SAGROW>
%         tblObjects(1) = getObjectsFromMask(0);

        % Make sure the sites with no crops found are labeled witn 0 num
        % crops
        missingObjectsList = Options.Objects(arrayfun(@(x) isempty(x.CropArea), Options.Objects));
        emptySites = unique([missingObjectsList.Site]);
        if ~isempty(emptySites)
            for e = emptySites
                Options.Crop(e).NumCrops = 0;
            end
        end
        
        for site = selectedSites
            siteIndex = getSiteIndex(file, site, metadata);
            well = getWellIndex(file, site, metadata);
            
            minCellSize = Options.Preprocessing(siteIndex).MinCellSize;
            if logical(Options.Crop(siteIndex).NumCrops)
                for cropIndex = 1:Options.Crop(siteIndex).NumCrops
                    objectIndex = getObjectIndex(file, siteIndex, cropIndex, Options.Objects);     

%                     fileKey = getFileKey(file, site, Options.Preprocessing(siteIndex).MaskChannel,...
%                         metadata, ExpSettings.Selected.Convention, well, cropIndex);
%                     stack = readFileToStack([ExpSettings.CroppedPath, fileKey, '.tif']);
%                     stack = double(stack);
%                     [~, mask, ~] = getCellMask(stack, Options.Preprocessing, siteIndex, cropIndex);
%                     
%                     fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);                    
%                     save([ExpSettings.MasksPath, fileKey, '-Mask.mat'], 'mask', '-v7.3');

                    fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                    uiprogressdlg(uifig, 'Title', 'Extracting Object Properties...', 'Message', fileKey, 'Value', index/length(indexes));
                    drawnow;

                    if ~isempty(Options.Objects(objectIndex).CropArea)

%                             Options.Objects(objectIndex).File = Options.Crop(siteIndex).File;
%                             Options.Objects(objectIndex).Site = Options.Crop(siteIndex).Site;
%                             Options.Objects(objectIndex).Crop = cropIndex;
%                             Options.Objects(objectIndex).FileKey = getFileKey(file, site,...
%                                Options.Preprocessing(siteIndex).MaskChannel, metadata,...
%                               ExpSettings.Selected.Convention, well, cropIndex);

                        fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                        load([ExpSettings.MasksPath, fileKey, '-Mask.mat'], 'mask');

                        fileKey = getFileKey(file, site, Options.Preprocessing(siteIndex).MaskChannel,...
                            metadata, ExpSettings.Selected.Convention, well, cropIndex);
                        stack = readFileToStack([ExpSettings.CroppedPath, fileKey, '.tiff']);

%                             Options.Objects(objectIndex).NumObjects2D = 1;
%                             Options.Objects(objectIndex).NumObjects3D = 1;

                        numObjects2D = zeros(size(mask,3),1);
                        for frame = 1:size(mask,3)
                            conComps2D = bwconncomp(mask(:,:,frame));
                            numObjects2D(frame,1) = conComps2D.NumObjects;                        
                        end
                        Options.Objects(objectIndex).NumObjects2D = max(numObjects2D, [], 'all');

                        conComps3D = bwconncomp(mask);
                        Options.Objects(objectIndex).NumObjects3D = conComps3D.NumObjects;

                        temp_objects = getObjectsFromMask(mask, minCellSize, stack);

                        if size(mask,3) == 1 % if single frame, save all together
                            if length(temp_objects) > 1 % More than one object in the current crop
                                % [~, indx] = min(abs(([temp_objects.BBX] -
                                % [temp_objects.BBY]))); % Keep only the
                                % one that is centered. Not good: keeps the
                                % one the most aligned with the diagonal.
                                [~, indx] = min(abs(Options.Objects(objectIndex).CropArea(3)/2 - [temp_objects.CentroidX]) + ...
                                    abs(Options.Objects(objectIndex).CropArea(4)/2 - [temp_objects.CentroidY])); % Keep only the one that is centered.
                                temp_objects = temp_objects(indx);
                                vectorizedMask = mask(:); vectorizedMask = false(size(vectorizedMask)); 
                                vectorizedMask(temp_objects.PixelIdxList{1,1}) = 1; 
                                mask = reshape(vectorizedMask, size(mask));

                                fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                                save([ExpSettings.MasksPath, fileKey, '-Mask.mat'], 'mask', '-v7.3');
       
                                % Saving outlined mask TIFF...
                                outlined_stack = cell(1,size(mask,3));
                                Stack2TIFF(uint16(mask(:,:,frame)), [ExpSettings.MasksPath, fileKey, '-Mask.tiff']);
                                for frame = 1:size(mask,3)       
                                    outlined_stack{frame} = DrawMaskOutline(stack(:,:,frame), mask(:,:,frame), [0,1,0]);
                                    imwrite(outlined_stack{frame}, ...
                                        [ExpSettings.MasksPath, fileKey,'-Outlined-Stack.tiff'], ...
                                        'WriteMode', WriteMode, 'Compression', 'none');                            
                                end                            
                                clear outlined_stack;
        
                                % Saving masked stack...
                                for channel = Options.Preprocessing(siteIndex).ChannelsToUse
                                    % Mask stack
                                    fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);            
                                    stack_path = [ExpSettings.CroppedPath, fileKey, '.tiff'];
                                    stack = readTIFFstack(stack_path);
                                    masked_stack = double(stack).*mask; masked_stack(masked_stack == 0) = NaN;          
                                    save([ExpSettings.MasksPath, fileKey, '-Masked-Stack.mat'], 'masked_stack', '-v7.3');
                                end                    
        
                                clear masked_stack;      


                            end

                            if temp_objects.NumObjects2D == 0
                                Options.Objects(objectIndex).NumObjects2D = 0;
                            end
                            
                            tblObjects(index) = temp_objects; % Append object
                            
                            Meta(index).MapIndex = indexes(index);
                            Meta(index).ExpID = ExpSettings.ExpID;
                            Meta(index).File = Options.Objects(objectIndex).File;
                            Meta(index).Site = Options.Objects(objectIndex).Site;
                            Meta(index).Crop = Options.Objects(objectIndex).Crop;
                            Meta(index).FileKey = Options.Objects(objectIndex).FileKey;
                            treatment = split(metadata(file).wellNames{well,1}, '-');
                            treatment = treatment(2,1); treatment = treatment{1,1};
                            Meta(index).Treatment = treatment;
                            % duration = split(metadata(file).wellNames{well,1}, '-');
                            % duration = duration(3,1); duration = duration{1,1};
                            % Meta(index).Duration = duration;
                            Meta(index).DateMod = datestr(now);
                            Options.Objects(objectIndex).DateModified = datestr(now);

                            if metadata(file).n_channels > 1
                                % Find number of nuclei using DAPI or HOECHST
                                Options.Preprocessing(siteIndex).NuclearChannel = 2;
                                fileKey = getFileKey(file, site, Options.Preprocessing(siteIndex).NuclearChannel,...
                                    metadata, ExpSettings.Selected.Convention, well);
                                nuclearStack = readFileToStack([ExpSettings.NoBackgroundPath, fileKey, '.tiff']);
                                nuclearStack = imcrop(nuclearStack, Options.Objects(objectIndex).CropArea);
                                nuclearOpts = Options.Preprocessing;
                                nuclearOpts(siteIndex).MaskChannel = Options.Preprocessing(siteIndex).NuclearChannel;
                                nuclearOpts(siteIndex).MinCellSize = 1000; % Minimum nucleus size
                                Options.Crop(siteIndex).MinNucleiSize = nuclearOpts(siteIndex).MinCellSize;
                                try
                                    [~, nuclearMask, ~] = getCellMask(nuclearStack, nuclearOpts, siteIndex, cropIndex);
                                    nuclearMask = mask.*nuclearMask;
                                    nuclearObjects = bwconncomp(nuclearMask);
                                    Options.Objects(objectIndex).NumNuclei = nuclearObjects.NumObjects;
                                catch
                                    Options.Objects(objectIndex).NumNuclei = NaN;
                                end
                            else
                                 Options.Objects(objectIndex).NumNuclei = NaN;
                            end

                            % Find boundary cells
                            [rows, cols] = ind2sub(size(mask), tblObjects(index).PixelIdxList{1,1});
                            rows = rows + floor(Options.Objects(objectIndex).CropArea(2));
                            cols = cols + floor(Options.Objects(objectIndex).CropArea(1));

                            if any(rows <= 2) || any(cols <= 2) ||...
                                    any(rows >= Options.Crop(siteIndex).OriginalSize(1)) || any(rows >= Options.Crop(siteIndex).OriginalSize(2)) || ...
                                    any(cols >= Options.Crop(siteIndex).OriginalSize(1)) || any(cols >= Options.Crop(siteIndex).OriginalSize(2))
                                Options.Objects(objectIndex).BoundaryCell = 1;
                            else
                                Options.Objects(objectIndex).BoundaryCell = 0;
                            end                                

                        else % If mask is a stack, save each stack's objects separately
                            
                                tblObjects = temp_objects; 
                                for frame = 1:length(tblObjects)
                                    Meta(frame).MapIndex = indexes(index);
                                    Meta(frame).File = Options.Objects(objectIndex).File;
                                    Meta(frame).Site = Options.Objects(objectIndex).Site;
                                    Meta(frame).Crop = Options.Objects(objectIndex).Crop;
                                    Meta(frame).FileKey = Options.Objects(objectIndex).FileKey;
%                                     treatment = split(metadata(file).wellNames{well,1}, '-');
%                                     treatment = treatment(2,1); treatment = treatment{1,1};
%                                     Meta(frame).Treatment = treatment;
%                                     duration = split(metadata(file).wellNames{well,1}, '-');
%                                     duration = duration(3,1); duration = duration{1,1};
%                                     Meta(frame).Duration = duration;
                                    Meta(frame).DateMod = datestr(now);
                                   
                                end
                                Options.Objects(objectIndex).DateModified = datestr(now);
                                Objects = table2struct([struct2table(Meta), struct2table(tblObjects)]);
                                fileKey = getFileKey(file, site, Options.Preprocessing(siteIndex).MaskChannel,...
                                    metadata, ExpSettings.Selected.Convention, well, cropIndex);
                                save([ExpSettings.ObjectsPath, fileKey, '-Objects.mat'], 'Objects', '-v7.3');    
                                
                        end


%                     else % If object is NOT labeled as retained
                        
%                         temp_objects = getObjectsFromMask(0);
%                         tblObjects(index) = temp_objects; % Append next object
%                         
%                         Options.Objects(objectIndex).NumObjects2D = NaN;
%                         Options.Objects(objectIndex).NumObjects3D = NaN;
%                         Options.Objects(objectIndex).NumNuclei = NaN;
%                         Options.Objects(objectIndex).BoundaryCell = NaN;
%                         
%                         Meta(index).MapIndex = indexes(index);
%                         Meta(index).File = Options.Objects(objectIndex).File;
%                         Meta(index).Site = Options.Objects(objectIndex).Site;
%                         Meta(index).Crop = Options.Objects(objectIndex).Crop;
%                         Meta(index).FileKey = Options.Objects(objectIndex).FileKey;
%                         Meta(index).DateMod = datestr(now);

                    end
                    index = index + 1;
                    % Meta = [];
                end
            else % If no crops for this site
                objectIndex = getObjectIndex(file, siteIndex, 1, Options.Objects);
                
                temp_objects = getObjectsFromMask(ones(Options.Crop(siteIndex).OriginalSize));
                tblObjects(index) = temp_objects; % Append next object              
                
                Options.Objects(objectIndex).NumObjects2D = NaN;
                Options.Objects(objectIndex).NumObjects3D = NaN;
                Options.Objects(objectIndex).NumNuclei = NaN;
                Options.Objects(objectIndex).BoundaryCell = NaN;

                Meta(index).MapIndex = indexes(index);
                Meta(index).ExpID = ExpSettings.ExpID;
                Meta(index).File = Options.Objects(objectIndex).File;
                Meta(index).Site = Options.Objects(objectIndex).Site;
                Meta(index).Crop = Options.Objects(objectIndex).Crop;
                Meta(index).FileKey = Options.Objects(objectIndex).FileKey;
                treatment = split(metadata(file).wellNames{well,1}, '-');
                treatment = treatment(2,1); treatment = treatment{1,1};
                Meta(index).Treatment = treatment;
                % duration = split(metadata(file).wellNames{well,1}, '-');
                % duration = duration(3,1); duration = duration{1,1};
                % Meta(index).Duration = duration;
                Meta(index).DateMod = datestr(now);                    

                index = index + 1;
            end           
        end
        
        if ~isempty(Meta)
            exportObjects = table2struct([struct2table(Meta), struct2table(tblObjects)]);
            exportObjects = rmfield(exportObjects, 'MapIndex');
        else
            exportObjects = tblObjects;
        end
        
        for i = 1:length(exportObjects)
            objectIndex = getObjectIndex(exportObjects(i).File, exportObjects(i).Site, exportObjects(i).Crop, Options.Objects);
            targetIndex = getObjectIndex(exportObjects(i).File, exportObjects(i).Site, exportObjects(i).Crop, Objects);
            if targetIndex > length(Objects)
                Objects(end+1) = exportObjects(i); %#ok<SAGROW>
            elseif length(objectIndex) == 1 && length(targetIndex) > 1
                while length(targetIndex) > 1
                    Objects(targetIndex(1)) = [];
                    targetIndex = getObjectIndex(exportObjects(i).File, exportObjects(i).Site, exportObjects(i).Crop, Objects);
                end                    
                Objects(targetIndex(1)) = exportObjects(i);
            else
                Objects(targetIndex) = exportObjects(i);
            end
        end              
        
        if ~isempty(fieldnames(Objects)) && length(Objects) > 1
            % Sort Objects based on File, then Site, then Crop
            sTbl = struct2table(Objects);
            [tbl, ~]  = sortrows(sTbl, {'File', 'Site', 'Crop'});
            Objects = table2struct(tbl); clear index tbl sTbl;
        end   
        
        
        if ~isfield(Objects, 'Treatment')
            fieldNames = fieldnames(Objects);
            newFieldNames = cell(size(fieldNames, 1) + 3, 1);
            for o = 1:length(Objects)
                Objects(o).ExpID = ExpSettings.ExpID;                
                well = getWellIndex(Objects(o).File, Objects(o).Site, metadata);
                Objects(o).Treatment = metadata(Objects(o).File).wellNames{well,1}(4:8);
                % Objects(o).Duration = metadata(Objects(o).File).wellNames{well,1}(12:13);
            end

            newFieldNames{1,1} = 'ExpID';
            newFieldNames{7,1} = 'Treatment';
            % newFieldNames{8,1} = 'Duration';
            f = 1;
            for nf = [2:6,9:length(newFieldNames)]
                newFieldNames{nf,1} = fieldNames{f,1};
                newFieldNames{nf,1} = fieldNames{f, 1};
                f = f + 1;
            end

            Objects = orderfields(Objects, newFieldNames);
        end
        
        fileKey = getFileKey(file, 0, 0, metadata, ExpSettings.Selected.Convention);
        uiprogressdlg(uifig, 'Title', 'Saving results...', 'Message', fileKey, 'Indeterminate', 'on');
        save([ExpSettings.ObjectsPath, fileKey, '-Objects.mat'], 'Objects', '-v7.3');
        objectsTable = struct2table(Objects);
        save([ExpSettings.ObjectsPath, fileKey, '-ObjectsTable.mat'], 'objectsTable', '-v7.3');
        
    end
    close(uifig);
    
    save([ExpSettings.CloudPath, ExpSettings.ExpID, '-ExpSettings.mat'], 'ExpSettings');
    save([ExpSettings.CloudPath, ExpSettings.ExpID, '-Options.mat'], 'Options', '-v7.3');    
            
end
