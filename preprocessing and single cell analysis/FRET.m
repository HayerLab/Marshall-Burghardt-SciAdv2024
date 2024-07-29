%% FRET Ratio
for f = 1 % Fake for loop to enable section collapse
    wb = waitbar(0, 'Starting FRET analysis...', 'Name', 'Eh! Fait FRET!');

    if exist(ExpSettings.CroppedPath, 'dir') && isfield(Options, 'Crop')
        sourcePath = ExpSettings.CroppedPath;
    else
        sourcePath = ExpSettings.NoBackgroundPath;
    end
    
    for file = ExpSettings.Selected.files_to_analyze
        bleach_raw_all = nan(metadata(file).n_frames,...
            length(ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))));    
        bleach_mean = nan(metadata(file).n_frames,...
            length(ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))));    
        corr = nan(metadata(file).n_frames,...
            length(ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))));    
        index = 1;
        for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
            % Select file
            siteIndex = getSiteIndex(file, site, metadata);            
            well = getWellIndex(file, site, metadata);
            
            if ~isfield(Options.FRET, 'refSites')
                for p = 1:length(Options.FRET)
                    Options.FRET(p).refSites = 1:metadata(Options.FRET(p).File).n_sites;
                end
            end
            
            if ~isfield(Options.FRET, 'FRETmode') || isempty(Options.FRET(siteIndex).FRETmode)
                Options.FRET(siteIndex).FRETmode = 'mean';
            end
            
            if strcmpi(sourcePath, ExpSettings.CroppedPath) && isfield(Options, 'Crop')
                nCrop = 1:Options.Crop(siteIndex).NumCrops;
            elseif ~isfield(Options, 'Crop')
                nCrop = 0;
            else
                nCrop = 1;
            end
            
            if strcmpi(Options.FRET(siteIndex).FilterType, 'Gaussian')
                h = fspecial(lower(Options.FRET(siteIndex).FilterType),...
                    Options.FRET(siteIndex).FilterSize,...
                    Options.FRET(siteIndex).FilterSize/6);
            elseif ~strcmpi(Options.FRET(siteIndex).FilterType, 'None')
                h = fspecial(lower(Options.FRET(siteIndex).FilterType),...
                    Options.FRET(siteIndex).FilterSize);
            end
            
            for cropIndex = nCrop
                
                fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                wb.Name = ['FRET - ', fileKey];
                waitbar(0, wb, 'Loading data...');
                
                fileKey = getFileKey(file, site, Options.FRET(siteIndex).DonorChannel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                if Options.FRET(siteIndex).MaskFRET
                    load([ExpSettings.MasksPath, fileKey, '-Masked-Stack.mat'], 'masked_stack');
                    DonorStack = masked_stack;                      
                    clear masked_stack;                
                else 
                    DonorStack = double(readTIFFstack([ExpSettings.NoBackgroundPath, fileKey, '.tif']));
                end
                

                fileKey = getFileKey(file, site, Options.FRET(siteIndex).AcceptorChannel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                if Options.FRET(siteIndex).MaskFRET            
                    load([ExpSettings.MasksPath, fileKey, '-Masked-Stack.mat'], 'masked_stack');
                    AcceptorStack = masked_stack;                
                    clear masked_stack;                
                else 
                    AcceptorStack = double(readTIFFstack([ExpSettings.NoBackgroundPath, fileKey, '.tif']));
                end
                
                fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                load([ExpSettings.MasksPath, fileKey, '-Mask.mat'], 'mask');

                imRatio = zeros(size(DonorStack));
                for frame = 1:size(DonorStack, 3)
                    if ~strcmpi(Options.FRET(siteIndex).FilterType, 'None')
%                         imCFP = imfilter(DonorStack(:,:,frame), h, 'replicate');
                        imCFP = ndnanfilter(DonorStack(:,:,frame), h, 'replicate');
                        imCFP(imCFP == 0) = NaN;

                        imYFP = ndnanfilter(AcceptorStack(:,:,frame), h, 'replicate');
                        imYFP(imYFP == 0) = NaN;
                    else
                        imCFP = DonorStack(:,:,frame); imCFP(imCFP == 0) = NaN;
                        imYFP = AcceptorStack(:,:,frame); imYFP(imYFP == 0) = NaN;
                    end                                

                    imRatio(:,:,frame) = imYFP./imCFP;

                    if Options.FRET(siteIndex).MaskFRET
                        imRatio(:,:,frame) = double(imRatio(:,:,frame)).*(mask(:,:,frame));
                    else
                        imRatio(:,:,frame) = double(imRatio(:,:,frame));
                    end

                    imRatio(imRatio == 0) = NaN;
                    tempRatio = imRatio(:,:,frame); tempMask = mask(:,:,frame); tempRatio(~tempMask) = NaN;
    %                 bleach_raw_all(frame, ExpSettings.Selected.sites(file,:) == site) = mean(tempRatio, "all", "omitnan");
%                     if sum(ismember(Options.FRET(siteIndex).Site, Options.FRET(siteIndex).refSites), 'all')
                    switch Options.FRET(siteIndex).FRETmode
                        case 'mean'
                            bleach_raw_all(frame, index) = mean(tempRatio, "all", "omitnan");
                        case 'median'
                            bleach_raw_all(frame, index) = median(tempRatio, "all", "omitnan");
                    end
%                     end
                    waitbar(frame/metadata(file).n_frames, wb, 'Calculating FRET...');
                end
%             bleach_raw_all(:,site) = bleach_raw_all(:,site)./median(bleach_raw_all(:, site), 1, 'omitnan');
            
                fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                save([ExpSettings.FRET, fileKey, '-FRET-Raw.mat'],'imRatio', '-v7.3');                

                index = index + 1;
            end
            
        end
        bleach_raw_all(bleach_raw_all == 0) = NaN;
        
        fileKey = getFileKey(file, 0, 0, metadata, ExpSettings.Selected.Convention);
        save([ExpSettings.FRET, fileKey, '-FRET-Raw.mat'],'bleach_raw_all', 'bleach_mean', '-v7.3');

        clear AcceptorStack DonorStack        
         
        %% Photobleach correction for FRET ratio files
        index = 1;
        if ~ishandle(wb)
            wb = waitbar(0, 'Starting FRET analysis...', 'Name', 'Eh! Fait FRET!');
        end
        
        for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
            % Select file
            siteIndex = getSiteIndex(file, site, metadata);
            well = getWellIndex(file, site, metadata);
            
%             switch Options.FRET(siteIndex).BleachCorrectPerFile                   
%                 case 1 % Gather all sites of each file to calculate the bleach curve.
%                     sitesRange = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file));
%                 case 0 % Fit bleaching model to each site.
%                     sitesRange = site;
%             end                  

            if ~isfield(Options.FRET, 'refFrames')
                for p = 1:length(Options.FRET)
                    Options.FRET(p).refFrames = 1:metadata(Options.FRET(p).File).n_frames;
                end
            end

            if strcmpi(sourcePath, ExpSettings.CroppedPath) && isfield(Options, 'Crop')
                nCrop = 1:Options.Crop(siteIndex).NumCrops;
            elseif ~isfield(Options, 'Crop')
                nCrop = 0;
            else
                nCrop = 1;
            end

            for cropIndex = nCrop                    
                fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);            
                wb.Name = ['FRET - ', fileKey];
                waitbar(0, wb, 'Fitting bleach model...');
                          
                
                index = index +1;
                
            end
%             [bleach_mean(:, index), fitpara, corr] = estimatePhotobleaching(bleach_raw_all(:, sitesRange), 'median',...
%                 Options.FRET(siteIndex).BleachingModel, 0);
        end
        
        if file > 1 && isfield(Options, 'Objects')
            [bleach_mean, fitpara, corr] = estimatePhotobleaching(bleach_raw_all(Options.FRET(siteIndex).refFrames,...
                getObjectIndex(file, Options.FRET(siteIndex).refSites(1), 1, Options.Objects) - getObjectIndex(file-1, metadata(file-1).n_sites,...
                    Options.Crop(getSiteIndex(file-1, Options.FRET(siteIndex).refSites(end), metadata)).NumCrops, Options.Objects):...
                getObjectIndex(file, Options.FRET(siteIndex).refSites(end),...
                    Options.Crop(getSiteIndex(file, Options.FRET(siteIndex).refSites(end), metadata)).NumCrops,...
                    Options.Objects) - getObjectIndex(file-1, metadata(file-1).n_sites,...
                    Options.Crop(getSiteIndex(file-1, Options.FRET(siteIndex).refSites(end), metadata)).NumCrops, Options.Objects)),...
                    'median', Options.FRET(siteIndex).BleachingModel, 1);
        elseif isfield(Options, 'Objects')
            [bleach_mean, fitpara, corr] = estimatePhotobleaching(bleach_raw_all(Options.FRET(siteIndex).refFrames,...
                getObjectIndex(file, Options.FRET(siteIndex).refSites(1), 1, Options.Objects):...
                getObjectIndex(file, Options.FRET(siteIndex).refSites(end),...
                    Options.Crop(getSiteIndex(file, Options.FRET(siteIndex).refSites(end), metadata)).NumCrops,...
                    Options.Objects)),...
                    'median', Options.FRET(siteIndex).BleachingModel, 1);
        else
            [bleach_mean, fitpara, corr] = estimatePhotobleaching(bleach_raw_all(Options.FRET(siteIndex).refFrames, Options.FRET(siteIndex).refSites),...
                    'median', Options.FRET(siteIndex).BleachingModel, 1);
        end
        


%         % Display results for single site. 
%         if Options.Preprocessing(siteIndex).PlotBleaching
%             [fig, ax] = plotBleachingFit(normalized_masked_mean, corr_norm);
%             hold(ax, 'on'); 
%             ylabel(ax, 'Normalized mean intensity of the masked area');
%             ttl = ax.Title; ttl.String = [ttl.String, ' for site ' fileKey]; ax.Title = ttl;
%             plt = plot(ax, 1:size(mean_bc_stack,1), (mean_bc_stack./median(mean_bc_stack, 1, "omitnan")), 'k');
%             plt.DisplayName = ['Bleach-corrected (', Options.Preprocessing(siteIndex).BleachingModel, ' model)'];
%             drawnow;
%             if Options.Preprocessing(siteIndex).SaveBleachingPlot
            fileKey = getFileKey(file, 0, 0, metadata, ExpSettings.Selected.Convention, well, 0);
            fig = gcf; ax = gca;
            savefig(fig, [ExpSettings.FRET, fileKey, '-Bleaching-Profiles.fig']);
%             exportgraphics(ax, [ExpSettings.FRET, fileKey, '-Bleaching-Profiles.png']);
%             end
%         end

        % If fit was done one a subset of frames:
        if size(corr, 1) < metadata(file).n_frames
            corr = repmat(median(corr, 1, 'omitnan'), [metadata(file).n_frames, 1]);
        end
        
        siteIndex = getSiteIndex(file, 1, metadata);
        if ~isfield(Options.FRET, 'NormalizeWRT')
            for p = 1:length(Options.FRET)
                Options.FRET(p).NormalizeWRT = 'Mean of refSites';
            end
        end
        switch Options.FRET(siteIndex).NormalizeWRT
            case 'Mean of refSites'
                corr_norm = bleach_mean;
            case 'Fitted trend'
                corr_norm = corr;
            case 'Normalized fitted trend'
                corr_norm = corr./median(corr);
        end
        
        index = 1;
        for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
            siteIndex = getSiteIndex(file, site, metadata);
            well = getWellIndex(file, site, metadata);
            if strcmpi(sourcePath, ExpSettings.CroppedPath) && isfield(Options, 'Crop')
                nCrop = 1:Options.Crop(siteIndex).NumCrops;
            elseif ~isfield(Options, 'Crop')
                nCrop = 0;
            else
                nCrop = 1;
            end         
            
            for cropIndex = nCrop
                fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                siteIndex = getSiteIndex(file, site, metadata);
                wb.Name = ['FRET - ', fileKey];
                waitbar(0, wb, 'Bleach correcting FRET...');
                load([ExpSettings.FRET, fileKey, '-FRET-Raw.mat'], 'imRatio');
                FRETdata = zeros(size(imRatio));
%                 meanFRET = zeros(size(imRatio,3), 1);
            
                for frame = 1:metadata(file).n_frames
                    FRETdata(:,:,frame) = double(imRatio(:,:,frame)./(corr_norm(frame, 1)));
                    switch Options.FRET(siteIndex).FRETmode
                        case 'mean'
                            trendFRET(frame, index) = mean(FRETdata(:,:,frame), "all", "omitnan"); %#ok<SAGROW>
                        case 'median'
                            trendFRET(frame, index) = median(FRETdata(:,:,frame), "all", "omitnan"); %#ok<SAGROW>
                    end
                    waitbar(frame/metadata(file).n_frames, wb);
                end
                
%                 [fig, ax] = plotBleachingFit(bleach_raw_all, corr);
%                 hold(ax, 'on'); 
%                 ylabel(ax, 'Normalized mean intensity of the masked area');                    
%                 ttl = ax.Title; ttl.String = [ttl.String, ' for site ' fileKey]; 
%                 plt = plot(ax, 1:size(meanFRET,1), (meanFRET), 'k');
%                 plt.DisplayName = ['Bleach-corrected (', Options.Preprocessing(siteIndex).BleachingModel, ' model)'];
                
                wb.Name = ['FRET - ', fileKey];
                waitbar(1, wb, 'Generating RGB & TIFF Stacks... ');
                if ~isfield(Options.FRET, 'colorRange') || isempty(Options.FRET(siteIndex).colorRange)
                    colorRange = [round(prctile(FRETdata, 0.01, "all"), 1), round(prctile(FRETdata, 99.99, "all"), 1)];
                    range = max(1-colorRange(1), colorRange(2)-1);
                    colorRange = [1-range, 1+range];
                    if sign(colorRange(1)) < 0; colorRange(1) = 0; end
                    Options.FRET(siteIndex).colorRange = colorRange;
                else
                    colorRange = Options.FRET(siteIndex).colorRange;
                end
                
                [~, colorRange] = Ratio2RGBTIFF(FRETdata,...
                    [ExpSettings.FRET, fileKey, '-FRET-RGB-',...
                    Options.FRET(siteIndex).FilterType(1),...
                    num2str(Options.FRET(siteIndex).FilterSize),...
                    '-', num2str(colorRange(1)), '-', num2str(colorRange(2)),'.tif'],...
                    Options.FRET(siteIndex).colorRange);
                
                FRET_TIFFStack = uint16(round(double(FRETdata.*10000)));
                FRET_TIFFStack(isnan(FRET_TIFFStack)) = 0;
                Stack2TIFF(FRET_TIFFStack,...
                    [ExpSettings.FRET, fileKey, '-FRET-TIFF-',...
                    Options.FRET(siteIndex).FilterType(1),...
                    num2str(Options.FRET(siteIndex).FilterSize),...
                    '-', num2str(colorRange(1)), '-', num2str(colorRange(2)),'.tif']);
                
                save([ExpSettings.FRET, fileKey, '-FRET.mat'],'FRETdata', 'colorRange', '-v7.3');
                Options.FRET(siteIndex).DateModified = datestr(now);
                waitbar(1, wb, 'Done!!'); 
                
                index = index + 1;
            end
        end
        fileKey = getFileKey(file, 0, 0, metadata, ExpSettings.Selected.Convention);
        save([ExpSettings.FRET, fileKey, '-NormalizedFRET.mat'], 'trendFRET', '-v7.3');
%         [fig, ax] = plotBleachingFit(bleach_raw_all, corr);
%         hold(ax, 'on');
%         ylabel(ax, 'Normalized mean intensity of the masked area');
%         ttl = ax.Title; ttl.String = [ttl.String, ' for site ' fileKey]; 
%         plt = plot(ax, 1:size(meanFRET,1), (meanFRET), 'k');
%         plt.DisplayName = ['Bleach-corrected (', Options.Preprocessing(siteIndex).BleachingModel, ' model)'];

        clear imRatio imCFP imYFP RGBstack bc_stack mean_bc_stack;
        
    end
    
    wb.Name = 'FRET';
    waitbar(1, wb, 'Saving options... ');
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-Options.mat'], 'Options');
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-ExpSettings.mat'], 'ExpSettings');

    
    %% Generate movies with custom color range: 
    % Change color Range in Options.FRET(siteIndex).colorRange first.
    
%     waitbar(0, wb, 'Generating RGB TIFF Stacks... ');  
%     index = 0;
%     for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
%         siteIndex = getSiteIndex(file, site, metadata);
%         well = getWellIndex(file, site, metadata);
%         if strcmpi(sourcePath, ExpSettings.CroppedPath) && isfield(Options, 'Crop')
%             nCrop = 1:Options.Crop(siteIndex).NumCrops;
%         elseif ~isfield(Options, 'Crop')
%             nCrop = 0;
%         else
%             nCrop = 1;
%         end
%         
%         for cropIndex = nCrop
%             fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
%             index = index +1; 
%             wb.Name = ['FRET - ', fileKey];
%             waitbar(index/length(ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))),...
%                 wb, ['Generating RGB TIFF Stacks for ', fileKey]);        
%             siteIndex = getSiteIndex(file, site, metadata);
%             load([ExpSettings.FRET, fileKey, '-FRET.mat'], 'FRETdata'); 
%             colorRange = Options.FRET(siteIndex).colorRange;
%             [~, colorRange] = Ratio2RGBTIFF(FRETdata,...
%                 [ExpSettings.FRET, fileKey, '-FRET-RGB-',...
%                     Options.FRET(siteIndex).FilterType(1),...
%                     num2str(Options.FRET(siteIndex).FilterSize),...
%                     '-', num2str(colorRange(1)), '-', num2str(colorRange(2)),'.tif'],...
%                 Options.FRET(siteIndex).colorRange);
%             
%                 FRET_TIFFStack = uint16(round(double(FRETdata.*10000)));
%                 FRET_TIFFStack(isnan(FRET_TIFFStack)) = 0;
%                 Stack2TIFF(FRET_TIFFStack,...
%                     [ExpSettings.FRET, fileKey, '-FRET-TIFF-',...
%                     Options.FRET(siteIndex).FilterType(1),...
%                     num2str(Options.FRET(siteIndex).FilterSize),...
%                     '-', num2str(colorRange(1)), '-', num2str(colorRange(2)),'.tif']);
%             
%             Options.FRET(siteIndex).DateModified = datestr(now);
%         end
%     end
%     wb.Name = 'FRET';
%     waitbar(1, wb, 'Saving options... ');      
%     save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-Options.mat'], 'Options');
%     save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-ExpSettings.mat'], 'ExpSettings');
    
    waitbar(1, wb,' Done!');  
    pause(3);
    close(wb);
end
