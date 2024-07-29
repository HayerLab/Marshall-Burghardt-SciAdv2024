%% Preprocessing

%% Check requirements
% Make sure that file folder paths and required variables are loaded.
if ~exist('metadata','var')
    disp('Warning: missing variable metadata. Please select it manually.');
    [metadata_filename, metadata_path] = uigetfile('*.mat');
    load(fullfile(metadata_path, metadata_filename), "-mat", 'metadata');
elseif ~exist('ExpSettings', 'var')
    ExpSettings = generateExpSettings(metadata);
end

%% Select files
if ~isfield(ExpSettings.Selected, 'sites')
    ExpSettings = selectWhich(metadata, ExpSettings, 'custom');
end
save([ExpSettings.CloudPath, ExpSettings.ExpID, '-ExpSettings']);

%% Choose Options
if ~exist('Options', 'var')
    Options = selectOptions(metadata, ExpSettings);
end
save([ExpSettings.CloudPath, ExpSettings.ExpID, '-Options.mat'], 'Options');
save([ExpSettings.CloudPath, ExpSettings.ExpID, '-ExpSettings.mat'], 'ExpSettings');

%% Convert ND2 to TIFF Stacks
nd2TIFFs(ExpSettings.RootPath, metadata, 0, 1, 'custom', ExpSettings);

%% Alignment
for a = 1 % Fake for loop to enable section collapse
% [optimizer, metric] = imregconfig('multimodal');
% alginmentStack = readTIFFstack('E:\Rodrigo\230518\Alignment\230518-20X-AVG-F1-EqStack.tif');
% tform = imregtform(YFP, CFP, 'rigid', optimizer, metric);
% figure; imshowpair(YFP, CFP);


    if isfield(Options, 'FRET') % Alignment is meant for FRET analysis        
        %% Generate average images across sites for channel alignment
        ExpSettings.AlignmentPath = [ExpSettings.RootPath, 'Alignment', filesep];
        if ~exist(ExpSettings.AlignmentPath, "dir")
            mkdir(ExpSettings.AlignmentPath);
        end
        
        frame = 1; % Choose a random frame
        index = 1;
        for channel = min(Options.FRET(1).AcceptorChannel, Options.FRET(1).DonorChannel) : max(Options.FRET(1).AcceptorChannel, Options.FRET(1).DonorChannel)%metadata(1).n_channels    
            for file = ExpSettings.Selected.files_to_analyze %length(metadata) % Alignment uses all files
                for site = 1:metadata(file).n_sites % and sites available for feature-rich images
                    if metadata(file).n_wells > 1
                        well = find(site <= cumsum(metadata(file).wells(:,3)), 1);
                    else
                        well = 1;
                    end
                    fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);

                    CH_Stack(:,:,index) = imread([ExpSettings.TIFFStacksPath, fileKey, '.tif'], frame); %#ok<SAGROW>
                    index = index + 1;
                end
            end
            CH_AVG = uint16(mean(CH_Stack,3));
            clear CH_Stack;

            fileKey = getFileKey(file, 0, channel, metadata, ExpSettings.Selected.Convention, 0);
            fileKey = [fileKey(1:6), fileKey(10:end)];
            imwrite(CH_AVG, [ExpSettings.AlignmentPath, fileKey, '-AVG-F', num2str(frame),'.tif'],'TIFF','Compression','None');
        end

        %% Generate alignment stack 
            alignStack = zeros(size(CH_AVG,1), size(CH_AVG,2), 2);
            nStack = zeros(size(CH_AVG,1), size(CH_AVG,2), 2); 
%             clear CH_AVG;

            for channel = min(Options.FRET(1).AcceptorChannel, Options.FRET(1).DonorChannel) : max(Options.FRET(1).AcceptorChannel, Options.FRET(1).DonorChannel)
                fileKey = getFileKey(file, 0, channel, metadata, ExpSettings.Selected.Convention, 0);
                fileKey = [fileKey(1:6), fileKey(10:end)];
                alignStack(:,:,channel) = imread([ExpSettings.AlignmentPath, fileKey, '-AVG-F', num2str(frame),'.tif']);
                nStack(:,:,channel) = normalize(double(imread([ExpSettings.AlignmentPath, fileKey, '-AVG-F', num2str(frame),'.tif'])));
            end

        %% Calculate channel alignment from average images 
        [pX, pY, dxMat, dyMat] = dualviewComputeAlignmentFromGridImages(alignStack);
        save([ExpSettings.AlignmentPath, ExpSettings.ExpID, ' - Alignment Parameters.mat'],'pX','pY','dxMat','dyMat', 'alignStack', 'frame');

        binning = 1; % relevant if alingment images and data images were acquired using distinct binning settings
        imaligned = dualviewAlignFromFittedSurface(alignStack,pX,pY,binning);
        
        %%

        figure;
        subplot(1,2,1); imagesc(dxMat); axis image; colorbar; title('Displacements in X');
        subplot(1,2,2); imagesc(dyMat); axis image; colorbar; title('Displacements in Y');
        sgtitle('Alignment displacement matrices');

        figure; ax = gca; 
        x = linspace(1, size(dxMat,2),size(dxMat,2));    x = repmat(x, [size(dxMat, 1), 1]);
        y = linspace(1, size(dxMat,1), size(dxMat,1))';    y = repmat(y, [1, size(dxMat, 2)]);
        try
            quiverc(ax, x, y, dxMat, dyMat, 1); axis ij; axis image;
%             quiverColor2D(x, y, dxMat, dyMat, 'Parent', ax); axis ij; axis image;
        catch
            quiver(ax, x, y, dxMat, dyMat, 1); axis ij; axis image;
        end
        ax.Color = 'k'; title('Alignment displacement quiver');

%         figure; 
%         subplot(1,2,1); imshowpair(nStack(:,:,1), nStack(:,:,2), 'diff'); title('Before alignment');
%         subplot(1,2,2); imshowpair(imaligned(:,:,1), imaligned(:,:,2), 'diff'); title('After alignment');
%         sgtitle('Normalized differences between channels');    

        figure;
        subplot(1,2,1), showImagesMergeChannels(alignStack(:,:,1),alignStack(:,:,2)); title('YFP/CFP merge before alignment');
        subplot(1,2,2), showImagesMergeChannels(imaligned(:,:,1),imaligned(:,:,2)); title('YFP/CFP merge after alignment');
        sgtitle('YFP/CFP Merged - Channel alignment using multi-sites image');
        
        clear alignStack

        %% Test alignment on sample data files
        
        for channel = 1:2
            if metadata(file).n_wells > 1
                well = find(site <= cumsum(metadata(file).wells(:,3)), 1);
            else
                well =1;
            end
            
            fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
    %         disp(fileKey)

            stack_path = [ExpSettings.TIFFStacksPath, fileKey, '.tif'];
            warning off
            stack = readTIFFstack(stack_path);

            imstack(:,:,channel) = double(stack(:,:,frame));%#ok<SAGROW>
        end

        imaligned = dualviewAlignFromFittedSurface(imstack(:,:,1:2), pX, pY, 1);
        imaligned = makeSquare(imaligned);

        figure;
        subplot(1,2,1), showImagesMergeChannels(imstack(:,:,1),imstack(:,:,2)); title('YFP/CFP merge before alignment');
        subplot(1,2,2), showImagesMergeChannels(imaligned(:,:,1),imaligned(:,:,2)); title('YFP/CFP merge after alignment');
        sgtitle('YFP/CFP Merged - Channel alignment using merged channels');
        
        clear imstack imaligned;

        %% Align data

        ExpSettings.AlignedTIFFStacksPath = [ExpSettings.TIFFStacksPath, 'Aligned', filesep];
        if ~exist('ExpSettings.AlignedTIFFStacksPath', 'dir')
            mkdir(ExpSettings.AlignedTIFFStacksPath)
        end

        for file = ExpSettings.Selected.files_to_analyze    
            imStack = zeros(metadata(file).x_pix, metadata(file).y_pix, 2);      

            for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))                              
                siteIndex = getSiteIndex(file, site, metadata);
                well = getWellIndex(file, site, metadata);
                
%                 imaligned = dualviewAlignFromFittedSurface(imStack{1}(:,:,1:2), pX, pY, 1);
%                 imCH1 = zeros(size(stack)); %comment out to
%                 "allow" the code to crop images to align them
%                 imCH2 = zeros(size(stack));
%                 imCH3 = zeros(size(stack));
%                 imCH4 = zeros(size(stack));

                
                for frame = 1:metadata(file).n_frames
                    channelsToAlign = Options.Preprocessing(siteIndex).ChannelsToUse; %sort([Options.FRET(siteIndex).DonorChannel, Options.FRET(siteIndex).AcceptorChannel]);
                    if frame == 1; WriteMode = "overwrite"; else; WriteMode = "append"; end
                    for channel = 1:metadata(file).n_channels
                        % Options.Preprocessing(siteIndex).ChannelsToUse % 1:2 % Assumes YFP = 1 & CFP = 2
                        fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);                                                
                        imStack(:,:,channel) = double(imread([ExpSettings.TIFFStacksPath, fileKey, '.tif'], frame));
                    end                                        

                    switch ceil(length(Options.Preprocessing(siteIndex).ChannelsToUse) / 2) %metadata(file).n_channels %
                        case 1
                            channelsToAlign = sort([Options.FRET(siteIndex).DonorChannel, Options.FRET(siteIndex).AcceptorChannel]);
                            
                            imaligned = dualviewAlignFromFittedSurface(imStack(:,:,channelsToAlign), pX, pY, 1);
                            imaligned = makeSquare(imaligned);
                            
                            imCH1 = imaligned(:,:,1);
                            fileKey = getFileKey(file, site, channelsToAlign(1), metadata, ExpSettings.Selected.Convention, well);                            
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 1), 'all')
                                imwrite(uint16(imCH1), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');
                            end

                            imCH2 = imaligned(:,:,2);  
                            fileKey = getFileKey(file, site, channelsToAlign(2), metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 2), 'all')
                                imwrite(uint16(imCH2), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');
                            end
                            
                        case 2
                            if metadata(file).n_channels ~=4; disp('Warning! Number of channels might be wrong!'); end
                            
                            imaligned = dualviewAlignFromFittedSurface(imStack(:,:,channelsToAlign), pX, pY, 1);
                            imaligned = makeSquare(imaligned);

                            imCH1 = imaligned(:,:,1);
                            fileKey = getFileKey(file, site, 1, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 1), 'all')
                                imwrite(uint16(imCH1), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');
                            end
                            
                            imCH2 = imaligned(:,:,2);                            
                            fileKey = getFileKey(file, site, 2, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 2), 'all')
                                imwrite(uint16(imCH2), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');            
                            end
                            
                            imaligned2 = dualviewAlignFromFittedSurface(imStack(:,:,3:4), pX, pY, 1);
                            imaligned2 = makeSquare(imaligned2);

                            imCH3 = imaligned2(:,:,1);
                            fileKey = getFileKey(file, site, 3, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 3), 'all')
                                imwrite(uint16(imCH3), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');
                            end

                            imCH4 = imaligned2(:,:,2);
                            fileKey = getFileKey(file, site, 4, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 4), 'all')
                                imwrite(uint16(imCH4), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');               
                            end
                        case 3
                            if metadata(file).n_channels ~= 6; disp('Warning! Number of channels might be wrong!'); end
                            
                            imaligned = dualviewAlignFromFittedSurface(imStack(:,:,channelsToAlign), pX, pY, 1);
                            imaligned = makeSquare(imaligned);

                            imCH1 = imaligned(:,:,1);
                            fileKey = getFileKey(file, site, 1, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 1), 'all')
                                imwrite(uint16(imCH1), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');
                            end
                            
                            imCH2 = imaligned(:,:,2);                            
                            fileKey = getFileKey(file, site, 2, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 2), 'all')
                                imwrite(uint16(imCH2), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');            
                            end
                            
                            imaligned2 = dualviewAlignFromFittedSurface(imStack(:,:,3:4), pX, pY, 1);
                            imaligned2 = makeSquare(imaligned2);

                            imCH3 = imaligned2(:,:,1);
                            fileKey = getFileKey(file, site, 3, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 3), 'all')
                                imwrite(uint16(imCH3), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');
                            end

                            imCH4 = imaligned2(:,:,2);
                            fileKey = getFileKey(file, site, 4, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 4), 'all')
                                imwrite(uint16(imCH4), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');               
                            end
                            
                            imaligned3 = dualviewAlignFromFittedSurface(imStack(:,:,5:6), pX, pY, 1);
                            imaligned3 = makeSquare(imaligned3);
                            
                            imCH3 = imaligned3(:,:,1);
                            fileKey = getFileKey(file, site, 3, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 5), 'all')
                                imwrite(uint16(imCH3), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');
                            end

                            imCH4 = imaligned3(:,:,2);
                            fileKey = getFileKey(file, site, 4, metadata, ExpSettings.Selected.Convention, well);
                            if sum(ismember(Options.Preprocessing(siteIndex).ChannelsToUse, 6), 'all')
                                imwrite(uint16(imCH4), [ExpSettings.AlignedTIFFStacksPath, fileKey, '.tif'],...
                                    'WriteMode', WriteMode,'Compression','none');               
                            end                            
                    end

                end                   

                disp(getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well));
                clear("imCH1", "imCH2", "imaligned", "imStack");
            end
        end
    else
        ExpSettings.AlignedTIFFStacksPath = [ExpSettings.TIFFStacksPath];
    end
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-ExpSettings.mat'], 'ExpSettings');
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-Options.mat'], 'Options');
end 

%% Jitter correction (Image registration)
for j = 1 % Fake for loop to enable section collapse
    ExpSettings.JitterPath = [ExpSettings.RootPath, 'Jitter', filesep];
    if ~exist(ExpSettings.JitterPath, "dir")
        mkdir(ExpSettings.JitterPath);
    end
    if requiredStep(Options.Preprocessing, 'JitterCorrection', ExpSettings.Selected, metadata) ||...
            requiredStep(Options.Preprocessing, 'StageCorrection', ExpSettings.Selected, metadata)
        %% Get Jitter from Stage Data
        if requiredStep(Options.Preprocessing, 'StageCorrection', ExpSettings.Selected, metadata)
            for file = ExpSettings.Selected.files_to_analyze
                for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
                    siteIndex = getSiteIndex(file, site, metadata);

                    % Define max jitter size
                    if ~isfield(Options.Preprocessing, 'MaxJitter')
                        jittersize = 20;
                    elseif Options.Preprocessing(siteIndex).MaxJitter > 1 %Max Jittersize specified in Options.
                        jittersize = Options.Preprocessing(siteIndex).MaxJitter; 
                    elseif Options.Preprocessing(siteIndex).MaxJitter < 1 % Max Jitter in percentage.
                        jittersize = round(metadata(file).x_pix * Options.Preprocessing(siteIndex).MaxJitter);
                    else
                        jittersize = 2;
                    end
                    
                    if ~isfield(Options.Preprocessing, 'JitterChannel')
                        for p = 1:length(Options.Preprocessing)
                            Options.Preprocessing(p).JitterChannel = 1;
                        end
%                         Options = changeOptions(Options, 'Preprocessing', 'JitterChannel', 1);
                    end
                    
                    Jitter = getJitterFromStageData(file, site, Options.Preprocessing(siteIndex).JitterChannel, metadata, ExpSettings.JitterPath);
                    disp(['Max jitter size = ', num2str(max(abs(min(Jitter, [], 1) - max(Jitter,[],1)), [], "all"))]);
                    if jittersize < max(abs(min(Jitter, [], 1) - max(Jitter,[],1)), [], "all")
                        disp(['Warning: Max jitter in site #', num2str(site), ' is ',...
                            num2str(max(abs(min(Jitter, [], 1) - max(Jitter,[],1)), [], "all")),...
                            ', which is bigger than ', num2str(jittersize),' pixels']);
                        jittersize = ceil(max(abs(min(Jitter, [], 1) - max(Jitter,[],1)), [], "all")/5)*5;
                        disp(['Adjusting jittersize to ', num2str(jittersize)]);                            
                    end
                    Options.Preprocessing(siteIndex).MaxJitter = jittersize;
                    
                end
            end
        end

        %% Apply Stage Jitter Correction
        if requiredStep(Options.Preprocessing, 'StageCorrection', ExpSettings.Selected, metadata)
            RegisteredStacks_path = [ExpSettings.AlignedTIFFStacksPath, 'Stage Registered', filesep];
            ExpSettings.RegisteredStacksPath = RegisteredStacks_path;
            if ~exist("RegisteredStacks_path", "dir")
                mkdir(RegisteredStacks_path);
            end

            disp('Correcting stage position jitter.');

            for file = ExpSettings.Selected.files_to_analyze
                for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
                    siteIndex = getSiteIndex(file, site, metadata);
                    
                    fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, 0);                    
                    [~, jitterFiles] = getFilenames(ExpSettings.JitterPath, [fileKey, '-StageJitter']);
                    if size(jitterFiles,1) > 1
                        disp(['Too many jitter files in ', ExpSettings.JitterPath, ' match "', fileKey,'"']);
                        break;
                    end
                    load([ExpSettings.JitterPath, jitterFiles{1}]);
                    jittersize = Options.Preprocessing(siteIndex).MaxJitter;
                    
                    if metadata(file).n_wells > 1
                        well = find(site <= cumsum(metadata(file).wells(:,3)), 1);
                    else
                        well =1;
                    end
                                        
                    for channel = Options.Preprocessing(siteIndex).ChannelsToUse

                        fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                        disp(fileKey)

                        stack_path = [ExpSettings.AlignedTIFFStacksPath, ...
                            fileKey, '.tif'];
                        warning off
                        stack = readTIFFstack(stack_path);                        
                        
                        RegisteredTiffStack = zeros(size(stack, 1) - 2*jittersize, size(stack, 2) - 2*jittersize, size(stack, 3), 'uint16');
                        for frame = 1:size(Jitter,1)
                            if frame == 1 || sum(Jitter(frame,:)) == 0
                                tempframecrop = CropJitter(stack(:,:,frame),jittersize,jittersize,jittersize,jittersize, 0, 0);
                            elseif sum(Jitter(frame,:)) ~= 0
                                tempframecrop = CropJitter(stack(:,:,frame),jittersize,jittersize,jittersize,jittersize, Jitter(frame,1), Jitter(frame,2));
                            end
                            RegisteredTiffStack(:,:,frame) = uint16(tempframecrop);
                            if frame == 1; WriteMode = "overwrite"; else; WriteMode = "append"; end
                            imwrite(RegisteredTiffStack(:,:,frame), [RegisteredStacks_path, fileKey, '.tif'],...
                                'WriteMode', WriteMode, 'Compression','none');
                        end
                    end
                    Options.Preprocessing(siteIndex).StageCorrection = datestr(now);
                    Options.Preprocessing(siteIndex).DateModified = datestr(now);
                end
            end
            clear RegisteredTiffStack;
            disp('Done!');
%         elseif isfield(Options, 'FRET')
%             ExpSettings.RegisteredStacksPath = ExpSettings.AlignedTIFFStacksPath;
        else
            ExpSettings.RegisteredStacksPath = ExpSettings.TIFFStacksPath;
        end

        %% Fourier-based Jitter Correction
        if requiredStep(Options.Preprocessing, 'JitterCorrection', ExpSettings.Selected, metadata)
            disp('Performing Jitter Correction...');

            if requiredStep(Options.Preprocessing, 'StageCorrection', ExpSettings.Selected, metadata)
                SourceStacks_path = [ExpSettings.AlignedTIFFStacksPath, 'Stage Registered', filesep];
                disp('...from Stage-Corrected TIFF Stacks');
            else
                SourceStacks_path = ExpSettings.AlignedTIFFStacksPath;            
                disp('...from original TIFF Stacks');
            end

            RegisteredStacks_path = [SourceStacks_path, 'Fourier Registered', filesep];
            ExpSettings.RegisteredStacksPath = RegisteredStacks_path;
            if ~exist("RegisteredStacks_path", "dir")
                mkdir(RegisteredStacks_path);
            end

            warning off
            for file = ExpSettings.Selected.files_to_analyze
                for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
                    siteIndex = getSiteIndex(file,site,metadata);
                    
                    if metadata(file).n_wells > 1
                        well = find(site <= cumsum(metadata(file).wells(:,3)), 1);
                    else
                        well =1;
                    end
                                        
                    for channel = [Options.Preprocessing(siteIndex).JitterChannel,...
                            setdiff(Options.Preprocessing(siteIndex).ChannelsToUse, Options.Preprocessing(siteIndex).JitterChannel)]
                        fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);                        

                        stack_path = [SourceStacks_path,...
                            fileKey, '.tif'];
                        clear stack;
                        stack = readTIFFstack(stack_path);

                        if channel == Options.Preprocessing(siteIndex).JitterChannel
                            fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well);
                            jittersize = 2;
                            [Jitter, newJittersize] = getJitterFromStack_FLB(stack_path, jittersize);
                            if jittersize ~= newJittersize
                                jittersize = newJittersize;
                                [Jitter, ~] = getJitterFromStack_FLB(stack_path, jittersize);
                            end
                            Options.Preprocessing(siteIndex).FourierJitter = jittersize;
                            fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, 0);
                            save([ExpSettings.JitterPath, fileKey, '-Jitter_FLB.mat'], 'Jitter');                

                            fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                            [RegisteredTiffStack, j] = RegisterStack_C(stack, Jitter, ...
                                [RegisteredStacks_path, fileKey], Options.Preprocessing(siteIndex).FourierJitter);
                            disp([fileKey, ' Max jitter: ', num2str(ceil(max(abs(min(Jitter, [], 1) - max(Jitter,[],1)), [], "all"))),...
                                ', cropped ', num2str(j), ' pixels.']);

                        else%if channel > Options.Preprocessing(siteIndex).JitterChannel
                            fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, 0);
                            load([ExpSettings.JitterPath, fileKey, '-Jitter_FLB.mat'], 'Jitter');

                            fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                            [RegisteredTiffStack, j] = RegisterStack_C(stack, Jitter,...
                                [RegisteredStacks_path, fileKey], Options.Preprocessing(siteIndex).FourierJitter);
                        end
                        disp(fileKey)                        
                    end
                    Options.Preprocessing(siteIndex).JitterCorrection = datestr(now);
                    Options.Preprocessing(siteIndex).DateModified = datestr(now);
                end
            end            
            disp('Done!');    
        end
    else
        ExpSettings.RegisteredStacksPath = ExpSettings.AlignedTIFFStacksPath;
    end
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-ExpSettings.mat'], 'ExpSettings');
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-Options.mat'], 'Options');
end

%% Background correction
for b = 1 % Fake for loop to enable section collapse
    if requiredStep(Options.Preprocessing, 'RemoveBackground', ExpSettings.Selected, metadata)
        %% Select background files and generate TIFF stacks        
        [~, files] = getFilenames([ExpSettings.BGPath, 'TIFF Stacks', filesep]);
        load([ExpSettings.BGPath, ExpSettings.ExpID,'-BGMetadata.mat'], 'BGmetadata');
        
        if isempty(files)
            disp(['There are no background TIFFS in ', ExpSettings.BGPath, 'TIFF Stacks', filesep, '!']);
            if length(BGmetadata) == 1
                bgmetadataFile = 1;                
            else
                open('BGmetadata');
                bgmetadataFile = input('Please enter the background file number to generate TIFFs');
            end
            nd2TIFFs(ExpSettings.BGPath, BGmetadata, 0, 1, 'custom', ExpSettings);
        else
            bgmetadataFile = 1;
            [~, files] = getFilenames([ExpSettings.BGPath, 'TIFF Stacks', filesep], 'AVG');
        end
        
        [~, bgmetadataFileName] = getFilenames(ExpSettings.BGPath, 'BGMetadata');
%         bgmetadataFileName = bgmetadataFileName(arrayfun(@(x) contains(x.name, 'BGMetadata'), bgmetadataFileName));
        if size(bgmetadataFileName, 2) == 1
            load([ExpSettings.BGPath, bgmetadataFileName{1}]);
%             load([ExpSettings.BGPath, bgmetadataFileName(1).name]);
        else
            [fileName, filePath] = uigetfile('*.mat');
            load([filePath, fileName]);
        end
        
        if isempty(files)
            %% Generate averaged background images
            index = 1;
%             for bgmetadataFile = 1:size(BGmetadata,1)
                bg_names = cell(BGmetadata(bgmetadataFile).n_channels,1);
                for channel = 1:BGmetadata(bgmetadataFile).n_channels
                    filenamecontains = split(char(BGmetadata(bgmetadataFile).channels{channel, 1}), ' ');
                    filenamecontains = filenamecontains{1,1};
                    if strcmpi(filenamecontains, 'NA'); continue; end
                    [~, bg_names{index,1}] = calculate_bg_img_rm_blobs([ExpSettings.BGPath, 'TIFF Stacks', filesep],...
                        filenamecontains);
                    index = index+1;
                end
%             end
        else
%             bgmetadataFile = 1;            
            bg_names = cell(BGmetadata(bgmetadataFile).n_channels,1);
            for channel = 1:metadata(ExpSettings.Selected.files_to_analyze(1)).n_channels
                if strcmpi( metadata(file).channels{channel,1}, 'NA'); continue; end
                bg_names{channel,1} = ['AVG-BG-', metadata(file).channels{channel,1}];
            end
        end
                
        if size(bg_names, 1) ~= metadata(ExpSettings.Selected.files_to_analyze(1)).n_channels
            disp("Warning: background files and data files don't have the same number of channels!");
        end
        
        for channel = 1:metadata(ExpSettings.Selected.files_to_analyze(1)).n_channels
            if isempty(bg_names{channel,1})
                continue
            else
                BG(:,:,channel) = double(imread([ExpSettings.BGPath, 'TIFF Stacks', filesep, bg_names{channel,1}, '.tif'])); %#ok<SAGROW>
            end
        end
                
        while size(BG, 3) < metadata(ExpSettings.Selected.files_to_analyze(1)).n_channels
            BG = cat(3, BG, zeros(size(BG(:,:,1))));
        end
        
        %% Align background images
%         if requiredStep(Options.Preprocessing, 'Alignment', ExpSettings.Selected, metadata)
        if  ~strcmpi(ExpSettings.AlignedTIFFStacksPath, ExpSettings.TIFFStacksPath)
            if ~exist('pX', 'var')
                load([ExpSettings.AlignmentPath, ExpSettings.ExpID, ' - Alignment Parameters.mat']);
            end
            binning = 1;

            if BGmetadata(bgmetadataFile).n_channels >= 2 % Assumes CH1 and  CH2 are for FRET on dual camera
                BG_Aligned = dualviewAlignFromFittedSurface(BG(:,:,1:2),pX,pY,binning);
                BG_Aligned_Square = makeSquare(BG_Aligned);

                CH1_AVG_BG_Aligned = BG_Aligned_Square(:,:,1);
                CH2_AVG_BG_Aligned = BG_Aligned_Square(:,:,2);
            end

            if BGmetadata(bgmetadataFile).n_channels >= 4 % Assumes dual camera for CH3 and CH4
                BG_Aligned = dualviewAlignFromFittedSurface(BG(:,:,3:4),pX,pY,binning);
                BG_Aligned_Square = makeSquare(BG_Aligned);

                CH3_AVG_BG_Aligned = BG_Aligned_Square(:,:,1);
                CH4_AVG_BG_Aligned = BG_Aligned_Square(:,:,2);
            end

            clear BG BG_Aligned BG_Aligned_Square;
        else
            switch metadata(file).n_channels
                case 1
                    CH1_AVG_BG_Aligned = BG;
                case 2
                    CH1_AVG_BG_Aligned = BG(:,:,1);
                    CH2_AVG_BG_Aligned = BG(:,:,2);
                case 3
                    CH1_AVG_BG_Aligned = BG(:,:,1);
                    CH2_AVG_BG_Aligned = BG(:,:,2);
                    CH3_AVG_BG_Aligned = BG(:,:,3);
                case 4
                    CH1_AVG_BG_Aligned = BG(:,:,1);
                    CH2_AVG_BG_Aligned = BG(:,:,2);
                    CH3_AVG_BG_Aligned = BG(:,:,3);
                    CH4_AVG_BG_Aligned = BG(:,:,4);                    
                otherwise
            end
        end


        %% Crop background images and subtract background
        ExpSettings.NoBackgroundPath = [ExpSettings.RegisteredStacksPath, 'Background corrected', filesep];
        if ~exist("ExpSettings.NoBackgroundPath", "dir")
            mkdir(ExpSettings.NoBackgroundPath);
        end

        for file = ExpSettings.Selected.files_to_analyze
            for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
                siteIndex = getSiteIndex(file, site, metadata);
                well = getWellIndex(file, site, metadata);
                for channel = Options.Preprocessing(siteIndex).ChannelsToUse
                    fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                    stack_path = [ExpSettings.RegisteredStacksPath,...
                        fileKey, '.tif'];
                    clear stack;
                    stack = readTIFFstack(stack_path);

                    if channel == 1
                        crop_size = (size(CH1_AVG_BG_Aligned, 1) - size(stack,1))/2;
                        
                        switch length(Options.Preprocessing(siteIndex).ChannelsToUse)%metadata(file).n_channels
                            case 1
                                CH1_AVG_BG_Aligned_cropped = cropBG(CH1_AVG_BG_Aligned, crop_size);
                                
                                AVG_BG_Aligned_cropped{1} = CH1_AVG_BG_Aligned_cropped; clear CH1_AVG_BG_Aligned_cropped;
                                
                            case 2
                                CH1_AVG_BG_Aligned_cropped = cropBG(CH1_AVG_BG_Aligned, crop_size);
                                CH2_AVG_BG_Aligned_cropped = cropBG(CH2_AVG_BG_Aligned, crop_size);

                                AVG_BG_Aligned_cropped{1} = CH1_AVG_BG_Aligned_cropped; clear CH1_AVG_BG_Aligned_cropped;
                                AVG_BG_Aligned_cropped{2} = CH2_AVG_BG_Aligned_cropped; clear CH2_AVG_BG_Aligned_cropped;
                                
                            case 3
                                CH1_AVG_BG_Aligned_cropped = cropBG(CH1_AVG_BG_Aligned, crop_size);
                                CH2_AVG_BG_Aligned_cropped = cropBG(CH2_AVG_BG_Aligned, crop_size);
                                CH3_AVG_BG_Aligned_cropped = cropBG(CH3_AVG_BG_Aligned, crop_size);

                                AVG_BG_Aligned_cropped{1} = CH1_AVG_BG_Aligned_cropped; clear CH1_AVG_BG_Aligned_cropped;
                                AVG_BG_Aligned_cropped{2} = CH2_AVG_BG_Aligned_cropped; clear CH2_AVG_BG_Aligned_cropped;
                                AVG_BG_Aligned_cropped{3} = CH3_AVG_BG_Aligned_cropped; clear CH3_AVG_BG_Aligned_cropped;
                                
                            case 4
                                CH1_AVG_BG_Aligned_cropped = cropBG(CH1_AVG_BG_Aligned, crop_size);
                                CH2_AVG_BG_Aligned_cropped = cropBG(CH2_AVG_BG_Aligned, crop_size);
                                CH3_AVG_BG_Aligned_cropped = cropBG(CH3_AVG_BG_Aligned, crop_size);
                                CH4_AVG_BG_Aligned_cropped = cropBG(CH4_AVG_BG_Aligned, crop_size);

                                AVG_BG_Aligned_cropped{1} = CH1_AVG_BG_Aligned_cropped; clear CH1_AVG_BG_Aligned_cropped;
                                AVG_BG_Aligned_cropped{2} = CH2_AVG_BG_Aligned_cropped; clear CH2_AVG_BG_Aligned_cropped;
                                AVG_BG_Aligned_cropped{3} = CH3_AVG_BG_Aligned_cropped; clear CH3_AVG_BG_Aligned_cropped;
                                AVG_BG_Aligned_cropped{4} = CH4_AVG_BG_Aligned_cropped; clear CH4_AVG_BG_Aligned_cropped
                        end

                        for frame = 1:metadata(file).n_frames
                            BG_mask(:,:,frame) = getBGMask(double(stack(:,:,frame))); %#ok<SAGROW>
                            imWO_BG = subBG(double(stack(:,:,frame)), ...
                                BG_mask(:,:,frame), AVG_BG_Aligned_cropped{channel});
                            if frame == 1; WriteMode = "overwrite"; else; WriteMode = "append"; end
                            imwrite(uint16(imWO_BG), [ExpSettings.NoBackgroundPath, fileKey, '.tif'],...
                                'WriteMode', WriteMode, 'Compression','none');                    
                        end
                    else
                        for frame = 1:metadata(file).n_frames                
                            imWO_BG = subBG(double(stack(:,:,frame)), ...
                                BG_mask(:,:,frame), AVG_BG_Aligned_cropped{channel});
                            if frame == 1; WriteMode = "overwrite"; else; WriteMode = "append"; end
                            imwrite(uint16(imWO_BG), [ExpSettings.NoBackgroundPath, fileKey, '.tif'],...
                                'WriteMode', WriteMode, 'Compression','none');
                        end                                
                    end               
                    fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                    disp(fileKey)                    
                end
                clear BG_mask imWO_BG AVG_BG_Aligned_cropped;
                siteIndex = getSiteIndex(file, site, metadata);
                Options.Preprocessing(siteIndex).RemoveBackground = datestr(now);
                Options.Preprocessing(siteIndex).DateModified = datestr(now);
            end
        end 
        disp('Background removal: Done!');
    else
        ExpSettings.NoBackgroundPath = ExpSettings.RegisteredStacksPath;
    end
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-ExpSettings.mat'], 'ExpSettings');
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-Options.mat'], 'Options', '-v7.3');
end

%% Cell cropping (optional)
for c = 1 % Fake for loop to enable section collapse

    if isfield(Options.Crop, 'Cropped')
    % Snippet to convert Options from the previous to the new architecture
    % (Options.Objects now goes crop by crop).
        for file = ExpSettings.Selected.files_to_analyze
            for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
                siteIndex = getSiteIndex(file, site, metadata);
                well = getWellIndex(file, site, metadata);
                cropList = getObjectIndex(file, site, 0, Options.Crop);
                for cropIndex = cropList
                    if ~isempty(Options.Crop(cropIndex).Cropped)
                        Options.Crop(cropIndex).NumCrops = size(Options.Crop(cropIndex).Cropped, 2);
                        for o = 1:size(Options.Crop(cropIndex).Cropped, 2)
                            objectIndex = getObjectIndex(file, site, o, Options.Objects);
                            Options.Objects(objectIndex).File = Options.Crop(cropIndex).File;
                            Options.Objects(objectIndex).Site = Options.Crop(cropIndex).Site;
                            Options.Objects(objectIndex).FileKey = getFileKey(file, site, 0, metadata,...
                                ExpSettings.Selected.Convention, well, o);
                            Options.Objects(objectIndex).CropArea = Options.Crop(cropIndex).Cropped(o).CropArea;
                            Options.Objects(objectIndex).Crop = o;
                        end
                        Options.Crop(cropIndex).CropMode = 'Manual';
                    end
                end
            end
        end
    end

    %% Decide whether to crop
    if ~isfield(Options.Preprocessing, 'Crop')
        cropSites = input('Do you want to crop? ');
        if cropSites
            for i = 1:length(Options.Preprocessing)
                Options.Preprocessing(i).Crop = 1;
            end
        else
            for i = 1:length(Options.Preprocessing)
                Options.Preprocessing(i).Crop = 0;
            end
        end
    end
    
    ExpSettings.CroppedPath = [ExpSettings.NoBackgroundPath, 'Cropped', filesep];
    if ~exist(ExpSettings.CroppedPath, 'dir'); mkdir(ExpSettings.CroppedPath); end

    if requiredStep(Options.Preprocessing, 'Crop', ExpSettings.Selected, metadata)
        if ~isfield(Options, 'Objects')
            Options.Objects = struct();
        end
        
        ExpSettings.ObjectsPath = [ExpSettings.RootPath, 'Objects', filesep];
        if ~exist(ExpSettings.ObjectsPath, 'dir'); mkdir(ExpSettings.ObjectsPath); end
        
        if ~isfield(Options, 'Crop')
            Options.Crop(length(Options.Preprocessing)) = struct();
        end
        
        uifig = uifigure(); uifig.WindowStyle = 'alwaysontop';
        uiprogressdlg(uifig,'Title',"Chop chop! Let's crop!", 'Message', 'Loading data...', 'Indeterminate', 'on');
        
        if isfield(Options, 'Objects') && ~isempty(fieldnames(Options.Objects)) && length(Options.Objects) > 1
            % Sort Options.Objects based on File, then Site, then Crop
            sTbl = struct2table(Options.Objects);
            [tbl, index]  = sortrows(sTbl, {'File', 'Site', 'Crop'});
            Options.Objects = table2struct(tbl); clear index tbl;
        end
        
        for file = ExpSettings.Selected.files_to_analyze
            selectedSites = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file));
            verboseSelection = uiconfirm(uifig, 'Do you want to see the results? ',...
                "Chop chop! Let's crop!", 'Options', {'Yes, show me! ', 'No, thanks.', 'Nevermind!'},...
                'DefaultOption', 2, 'CancelOption', 3);
            verboseMode = strcmpi(verboseSelection, 'Yes, show me! ');
            for site = selectedSites
                cropSite = 0; cropIndex = 1;
                siteIndex = getSiteIndex(file, site, metadata);
                well = getWellIndex(file, site, metadata);
                objectIndex = getObjectIndex(file, siteIndex, cropIndex, Options.Objects);

                Options.Crop(siteIndex).File = Options.Preprocessing(siteIndex).File;
                Options.Crop(siteIndex).Site = Options.Preprocessing(siteIndex).Site;
                Options.Crop(siteIndex).FileKey = getFileKey(file, site, 0, metadata,...
                    ExpSettings.Selected.Convention, well);
                
                channel = Options.Preprocessing(siteIndex).MaskChannel;
                fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);                
                stack = readFileToStack([ExpSettings.NoBackgroundPath, fileKey, '.tiff']);

                % StackViewer(stack);
                
                fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well);
                if ~isfield(Options.Crop, 'CropMode') || site == selectedSites(1)
                     cropModeSelection = uiconfirm(uifig, ['Select Crop Mode for file:', fileKey],...
                        "Chop chop! Let's crop!", 'Options', {'Manual', 'Automatic', 'Polygon', 'Nevermind!'},...
                        'DefaultOption', 1, 'CancelOption', 4);
                     fileKey = getFileKey(file, 0, 0, metadata, ExpSettings.Selected.Convention, well);
                     useSame = uiconfirm(uifig, ['Use the same Crop Mode for all sites in file:', fileKey, '? '],...
                        "Chop chop! Let's crop!", 'Options', {'Yes!', 'No, let me choose for each.', 'Nevermind!'},...
                        'DefaultOption', 1, 'CancelOption', 3);
                     allSame = 1; globalSelection = cropModeSelection;
                elseif isfield(Options.Crop, 'NumCrops') && ~isempty(Options.Crop(siteIndex).NumCrops) && Options.Crop(siteIndex).NumCrops == 0
                    if allSame
                        cropModeSelection = globalSelection;
                    else
                        cropModeSelection = uiconfirm(uifig, ['Select Crop Mode for file:', fileKey],...
                            "Chop chop! Let's crop!", 'Options', {'Manual', 'Automatic', 'Polygon', 'Nevermind!'},...
                            'DefaultOption', 1, 'CancelOption', 4);
                        allSame = 0; globalSelection = [];
                    end
                elseif isfield(Options.Crop, 'CropMode') && isempty(Options.Crop(siteIndex).CropMode)
                    if allSame
                        cropModeSelection = globalSelection;
                    else                    
                         cropModeSelection = uiconfirm(uifig, ['Select Crop Mode for file:', fileKey],...
                            "Chop chop! Let's crop!", 'Options', {'Manual', 'Automatic', 'Polygon', 'Nevermind!'},...
                            'DefaultOption', 1, 'CancelOption', 4);
                        allSame = 0; globalSelection = [];
                    end
                else
                    cropModeSelection = Options.Crop(siteIndex).CropMode;
                end
                
                if  (isfield(Options.Crop, 'CropMode') && isempty(Options.Crop(siteIndex).CropMode)) || allSame
                    Options.Crop(siteIndex).CropMode = globalSelection;
                else
                    Options.Crop(siteIndex).CropMode = cropModeSelection;
                end
                clear cropModeSelection;
                
                Options.Objects(objectIndex).File = Options.Crop(siteIndex).File;
                Options.Objects(objectIndex).Site = Options.Crop(siteIndex).Site;
                Options.Objects(objectIndex).Crop = cropIndex;
                Options.Objects(objectIndex).FileKey = getFileKey(file, site, channel, metadata,...
                    ExpSettings.Selected.Convention, well, cropIndex);
                
                Options.Crop(siteIndex).OriginalSize = size(stack);
                Options.Crop(siteIndex).NumCrops = 0;

                switch Options.Crop(siteIndex).CropMode
                    case 'Manual'
                        fig = figure; ax = gca; hold(ax, 'on'); fig.WindowState = 'maximized';
%                         fig.Color = 'k'; 
%                         for frame = 1:size(stack, 3)
%                             imagesc(ax, stack(:,:,frame), [0, 300]); colormap('gray'); hold on;
%                             colorbar; title([fileKey, ' Frame # ', num2str(frame)]);
%                             axis image; drawnow;
%                         end
                        imagesc(ax, mean(stack, 3)); axis image; axis ij; colormap('hot');
                        title(ax, [fileKey, ' Mean intensity']);
                        maxShow = 300;
                        clim([0,maxShow]);
                        StackViewer(stack);

                        cropSelection = uiconfirm(uifig, 'Do you want to crop the image?',...
                            fileKey, 'Options', {'Yes, crop!', 'Nah thanks, this looks good.'}, 'DefaultOption', 1);
                        
                        objectIndex = getObjectIndex(file, siteIndex, 1, Options.Objects);
                        if strcmpi(cropSelection, 'Yes, crop!')
                            cropSite = 1; cropIndex = 1;
                            % Options.Objects(objectIndex).NumObjects2D = 0;
                        elseif strcmpi(cropSelection, 'Nah thanks, this looks good.') && isempty(Options.Objects(objectIndex).CropArea)
                            cropSite = 0;
                            Options.Objects(objectIndex).Retained = 0;
                            Options.Objects(objectIndex).DateModified = datestr(now);
                            Options.Objects(objectIndex).NumObjects2D = NaN;
                            Options.Objects(objectIndex).NumObjects3D = NaN;
                            Options.Objects(objectIndex).NumNuclei = NaN;
                            Options.Objects(objectIndex).BoundaryCell = NaN;                            
                        end
                        
                        if verboseMode && ishandle(fig)
                            imagesc(ax, mean(stack, 3)); axis image; axis ij; colorbar; colormap('gray');
                            title(ax, fileKey);
                            subtitle(ax, 'Manual cropping');
%                         ax.CLim(2) = ax.CLim(2)*0.4; caxis(ax.CLim);
                            clim([0,maxShow]);
                        end

                        while cropSite == 1
                            channel = Options.Preprocessing(siteIndex).MaskChannel;
                            fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                            refImage = mean(stack, 3); refImage(refImage > maxShow) = NaN;
                            [~, cropArea] = serimcropold(stack, refImage);
                            objectIndex = getObjectIndex(file, siteIndex, cropIndex, Options.Objects);
                            if objectIndex > length(Options.Objects)
                                Options.Objects(objectIndex) = Options.Objects(objectIndex-1);
                                Options.Objects(objectIndex).File = Options.Crop(siteIndex).File;
                                Options.Objects(objectIndex).Site = Options.Crop(siteIndex).Site;
                                Options.Objects(objectIndex).Crop = cropIndex;
                                Options.Objects(objectIndex).FileKey = getFileKey(file, site, channel, metadata,...
                                    ExpSettings.Selected.Convention, well, cropIndex);                                
                            end                            
                            Options.Objects(objectIndex).CropArea = cropArea;
                            Options.Objects(objectIndex).Retained = 1;
                            Options.Objects(objectIndex).DateModified = datestr(now);

                            cropSelection = uiconfirm(uifig, 'Do you want to crop another region in the same image?',...
                                fileKey, 'Options', {'Yes, crop again!', char("That's enough thanks!")}, 'DefaultOption', 1);
                            if strcmpi(cropSelection, 'Yes, crop again!')
                                cropSite = 1;
                                cropIndex = cropIndex + 1;
                            else
                                cropSite = 0;
                            end
                            
                            if verboseMode
                                rectangle('Position', Options.Objects(objectIndex).CropArea, 'EdgeColor', 'c', 'LineWidth', 1, 'Parent', ax);
                                text(Options.Objects(objectIndex).CropArea(1) + 15,... %+0.5*Options.Objects(objectIndex).CropArea(3),...
                                    Options.Objects(objectIndex).CropArea(2) + 20,... %0.5*Options.Objects(objectIndex).CropArea(4),...
                                    num2str(Options.Objects(objectIndex).Crop), 'Color', 'c', 'Parent', ax);
                            end
                            
                        end
                        Options.Crop(siteIndex).NumCrops = cropIndex;
                        Options.Preprocessing(siteIndex).Threshold = zeros(metadata(file).n_frames, Options.Crop(siteIndex).NumCrops);
                        
                        if verboseMode && ishandle(fig)
                            fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                            drawnow;
                            if ~exist([ExpSettings.TIFFStacksPath, 'Crop Schemes', filesep], 'dir'); mkdir([ExpSettings.TIFFStacksPath, 'Crop Schemes', filesep]); end
                            exportgraphics(fig, [ExpSettings.TIFFStacksPath, 'Crop Schemes', filesep, fileKey, '-CropScheme.tiff']);
                        end                         
                        
                        close(fig);
                    
                    case 'Automatic'
                        fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                        uiprogressdlg(uifig, 'Message', ['File ', fileKey, ' -> '], 'Value',...
                            (find(ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file)) == site)/...
                            length(ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))))); 
                        uifig.Name = 'Rough segmentation, object detection & cropping... ';

                        fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                        if verboseMode
                            fig = figure; ax = gca; hold(ax, 'on'); fig.WindowState = 'maximized';
                            imagesc(ax, mean(stack, 3)); axis image; axis ij; colorbar; colormap('gray');
                            title(ax, fileKey);
                            clim([50,500]);
                        end                        
                        
                        if strcmpi(Options.Crop(siteIndex).CropMode, 'Manual') && ~isempty(Options.Objects(objectIndex).CropArea)
                            Options.Crop(siteIndex).CropMethod = 'Manual';
                            if verboseMode
                                subtitle(ax, 'Manual cropping');
                                for cropIndex = 1:Options.Crop(siteIndex).NumCrops
                                    objectIndex = getObjectIndex(file, siteIndex, cropIndex, Options.Objects);
                                    rectangle('Position', Options.Objects(objectIndex).CropArea, 'EdgeColor', 'c', 'LineWidth', 1, 'Parent', ax);
                                    text(Options.Objects(objectIndex).CropArea(1) + 15,... +0.5*Options.Objects(objectIndex).CropArea(3),...
                                        Options.Objects(objectIndex).CropArea(2) + 20,... +0.5*Options.Objects(objectIndex).CropArea(4),...
                                        num2str(Options.Objects(objectIndex).Crop), 'Color', 'c', 'Parent', ax);
                                end
                            end
                            
                        else
                            if verboseMode
                                subtitle(ax, 'Automatic cropping');
                            end
                            MarginPx = 128;
                            Options.Crop(siteIndex).MarginPx = MarginPx;
                            AutoStruct = Options.Preprocessing;
                            if size(stack,3) >1
                                segImg = round(sum(imgaussfilt3(stack), 3, 'omitnan'));
                            else
                                segImg = round(imgaussfilt(stack, 3));
                            end
                            [~, mask, AutoStruct] = getCellMask(segImg, AutoStruct, siteIndex);
                            Options.Crop(siteIndex).CropMethod = AutoStruct(siteIndex).MaskMethod;
                            Options.Crop(siteIndex).MinThreshold = AutoStruct(siteIndex).minThreshold;
                            Options.Crop(siteIndex).MinCellSize = AutoStruct(siteIndex).MinCellSize;
                            objectsFound = getObjectsFromMask(mask, Options.Preprocessing(siteIndex).MinCellSize-1, stack);

                            for cropIndex = 1:length(objectsFound)
                                Options.Objects(objectIndex).File = Options.Crop(siteIndex).File;
                                Options.Objects(objectIndex).Site = Options.Crop(siteIndex).Site;
                                Options.Objects(objectIndex).Crop = cropIndex;
                                Options.Objects(objectIndex).FileKey = getFileKey(file, site, channel, metadata,...
                                    ExpSettings.Selected.Convention, well, cropIndex);                            
                                
                                Options.Objects(objectIndex).CropArea = [max(objectsFound(cropIndex).BBX - MarginPx, 1),...
                                    max(objectsFound(cropIndex).BBY - MarginPx, 1),...
                                    min((objectsFound(cropIndex).BBW + 2*MarginPx), (size(mask, 2) - objectsFound(cropIndex).BBX + MarginPx)),...
                                    min((objectsFound(cropIndex).BBH + 2*MarginPx), (size(mask, 1) - objectsFound(cropIndex).BBY + MarginPx))];
                                Options.Objects(objectIndex).DateModified = datestr(now);
    
                                if isempty(objectsFound) || (length(objectsFound) == 1 &&...
                                        ((objectsFound.BBX == 0.5 || objectsFound.BBY == 0.5) &&...
                                        (objectsFound.BBW == size(mask,2) || objectsFound.BBH == size(mask,1))))
                                    Options.Objects(objectIndex).CropArea = [];
                                    Options.Objects(objectIndex).Retained = 0;
                                    objectIndex = objectIndex + 1;
                                    continue
                                else
                                    Options.Objects(objectIndex).Retained = 1;
                                    if verboseMode
                                        rectangle('Position', Options.Objects(objectIndex).CropArea, 'EdgeColor', 'r', 'LineWidth', 1);
                                        text(Options.Objects(objectIndex).CropArea(1) + 15,... +0.5*Options.Objects(objectIndex).CropArea(3),...
                                            Options.Objects(objectIndex).CropArea(2) + 20,... +0.5*Options.Objects(objectIndex).CropArea(4),...
                                            num2str(Options.Objects(objectIndex).Crop), 'Color', 'm', 'Parent', ax);
                                    end
                                end
                                
                                objectIndex = objectIndex + 1;
                            end  
                            Options.Crop(siteIndex).NumCrops = cropIndex*Options.Objects(objectIndex-1).Retained;
                        end

                        uiprogressdlg(uifig, 'Message', ['File ', fileKey, ' -> ', num2str(cropIndex), ' crops.'], 'Indeterminate', 'on');
                        
                        if verboseMode
                            fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                            drawnow;
                            if ~exist([ExpSettings.TIFFStacksPath, 'Crop Schemes', filesep], 'dir'); mkdir([ExpSettings.TIFFStacksPath, 'Crop Schemes', filesep]); end
                            exportgraphics(fig, [ExpSettings.TIFFStacksPath, 'Crop Schemes', filesep, fileKey, '-CropScheme.tiff']);
                            close(fig);
                        end
                                          
                    case 'Polygon'
                        channel = Options.Preprocessing(siteIndex).MaskChannel;
                        fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);                          
                        fig = figure; ax = gca; hold(ax, 'on'); fig.WindowState = 'maximized';
%                         fig.Color = 'k'; 
%                         for frame = 1:size(stack, 3)
%                             imagesc(ax, stack(:,:,frame)); 
% %                             colormap('gray'); 
%                             hold on;
% %                             caxis([ax.CLim(1), ax.CLim(2)*0.2]);
%                             colorbar; title([fileKey, ' Frame # ', num2str(frame)]);
%                             axis image; axis ij; drawnow;
%                         end
                        imagesc(ax, mean(stack, 3)); axis image; axis ij; colormap('hot');
                        % imagesc(ax, stack(:,:,1)); axis image; axis ij; colormap('hot');
                        title(ax, [fileKey, ' Mean intensity']);
            %             StackViewer(stack);
                        maxShow = 300;
                        clim([0,maxShow]);

                        fig.WindowState = "maximized";

                        cropSelection = uiconfirm(uifig, 'Do you want to crop the image?',...
                            fileKey, 'Options', {'Yes, crop!', 'Nah thanks, this looks good.'}, 'DefaultOption', 1);
                        
                        objectIndex = getObjectIndex(file, siteIndex, 1, Options.Objects);
                        if strcmpi(cropSelection, 'Yes, crop!')
                            cropSite = 1; cropIndex = 1;
                            % Options.Objects(objectIndex).NumObjects2D = 0;
                        elseif strcmpi(cropSelection, 'Nah thanks, this looks good.') && isempty(Options.Objects(objectIndex).CropArea)
                            cropSite = 0;
                            Options.Objects(objectIndex).Retained = 0;
                            Options.Objects(objectIndex).DateModified = datestr(now);
                            Options.Objects(objectIndex).NumObjects2D = NaN;
                            Options.Objects(objectIndex).NumObjects3D = NaN;
                            Options.Objects(objectIndex).NumNuclei = NaN;
                            Options.Objects(objectIndex).BoundaryCell = NaN;                            
                        end
                        
                        if verboseMode && ishandle(fig)
                            imagesc(ax, mean(stack, 3)); axis image; axis ij; colorbar; colormap('gray');
                            title(ax, fileKey);
                            subtitle(ax, 'Polygon cropping');
%                         ax.CLim(2) = ax.CLim(2)*0.4; caxis(ax.CLim);
                            clim([0,maxShow]);
                        end

                        while cropSite == 1
                            channel = Options.Preprocessing(siteIndex).MaskChannel;
                            fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                            refImage = mean(stack, 3); refImage(refImage > maxShow) = NaN;
                            [cellMask, polygonVertices] = selectPolygonMask(refImage);
                            % cellMask = selectCellMaskFromAverageImage(sum(stack, 3));
                            % cellMask = selectCellMaskFromAverageImage(stack(:,:,1));
                            cellBBox = regionprops(cellMask, 'BoundingBox');
                            cropArea = cellBBox.BoundingBox;
                            
                            objectIndex = getObjectIndex(file, siteIndex, cropIndex, Options.Objects);
                            if objectIndex > length(Options.Objects)
                                Options.Objects(objectIndex) = Options.Objects(objectIndex-1);
                                Options.Objects(objectIndex).File = Options.Crop(siteIndex).File;
                                Options.Objects(objectIndex).Site = Options.Crop(siteIndex).Site;
                                Options.Objects(objectIndex).Crop = cropIndex;
                                Options.Objects(objectIndex).FileKey = getFileKey(file, site, channel, metadata,...
                                    ExpSettings.Selected.Convention, well, cropIndex);                                
                            end                            
                            Options.Objects(objectIndex).CropArea = cropArea;
                            Options.Objects(objectIndex).CropVertices = polygonVertices;
                            Options.Objects(objectIndex).Retained = 1;
                            Options.Objects(objectIndex).DateModified = datestr(now);

                            cropSelection = uiconfirm(uifig, 'Do you want to crop another region in the same image?',...
                                fileKey, 'Options', {'Yes, crop again!', char("That's enough thanks!")}, 'DefaultOption', 1);
                            if strcmpi(cropSelection, 'Yes, crop again!')
                                cropSite = 1;
                                cropIndex = cropIndex + 1;
                            else
                                cropSite = 0;
                            end
                            
                            if verboseMode
                                rectangle('Position', Options.Objects(objectIndex).CropArea, 'EdgeColor', 'c', 'LineWidth', 1.5, 'Parent', ax);
                                plot(polygonVertices.x, polygonVertices.y, '-c', 'LineWidth', 1, 'Parent', ax);
                                text(Options.Objects(objectIndex).CropArea(1) + 15,... %+0.5*Options.Objects(objectIndex).CropArea(3),...
                                    Options.Objects(objectIndex).CropArea(2) + 20,... %0.5*Options.Objects(objectIndex).CropArea(4),...
                                    num2str(Options.Objects(objectIndex).Crop), 'Color', 'c', 'Parent', ax);
                            end
                        end
                        Options.Crop(siteIndex).NumCrops = cropIndex;
                        Options.Preprocessing(siteIndex).Threshold = zeros(metadata(file).n_frames, Options.Crop(siteIndex).NumCrops);
                        
                        if verboseMode && ishandle(fig)
                            fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                            drawnow;
                            if ~exist([ExpSettings.TIFFStacksPath, 'Crop Schemes', filesep], 'dir'); mkdir([ExpSettings.TIFFStacksPath, 'Crop Schemes', filesep]); end
                            exportgraphics(fig, [ExpSettings.TIFFStacksPath, 'Crop Schemes', filesep, fileKey, '-CropScheme.tiff']);
                        end                         
                        
                        close(fig);                        
                        
                    case 'Nevermind!'
                        close(fig);
                        continue
                end
                clear Objects stack;       
            end
            
            failedCrops = Options.Objects(arrayfun(@(x) isempty(x.CropArea), Options.Objects));
            disp(['There are currently ', num2str(length(failedCrops)),...
                ' unsuccessful crops out of ', num2str(length(Options.Objects)),...
                ' (', num2str(round(100*length(failedCrops)/length(Options.Objects))),'% failed)']);
            failedSites = Options.Objects(find(arrayfun(@(x) isempty(x.CropArea), Options.Objects)));
            failedSites = [failedSites.Site];
            disp('Use the following code to select sites to fix: ');
            disp('ExpSettings.Selected.sites = zeros(size(ExpSettings.Selected.sites)); ExpSettings.Selected.sites(1,1:length(failedSites)) = failedSites; ExpSettings.Selected.n_sites_per_file = length(failedSites);');
            
            disp('Use the following code to select the recently cropped sites: ');
            disp('unique([Options.Objects(find(isnan([Options.Objects.NumObjects2D]))).Site]);');
            
            %% Execute crop action using selected Crop Areas
            if ~ishandle(uifig)
                uifig = uifigure(); uifig.WindowStyle = 'alwaysontop';                
            end
            uiprogressdlg(uifig,'Title',"Chop chop! Let's crop!", 'Message', 'Loading data...', 'Indeterminate', 'on');
            
            if isfield(Options, 'Objects') && ~isempty(fieldnames(Options.Objects)) && length(Options.Objects) > 1
                % Sort Options.Objects based on File, then Site, then Crop                
                sTbl = struct2table(Options.Objects);
                % isempty(sTbl(:, "File"))
                [tbl, index]  = sortrows(sTbl, {'File', 'Site', 'Crop'});
                Options.Objects = table2struct(tbl); clear index tbl;
            end
            warning  off

            remainingCrops = [];
            for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
                siteIndex = getSiteIndex(file, site, metadata);
                well = getWellIndex(file, site, metadata);

                if isfield(Options, 'Crop')
                    if ~exist('remainingCrops', 'var') || isempty(remainingCrops)
                        remainingCrops = 1:Options.Crop(siteIndex).NumCrops;
                    end
                else
                    remainingCrops = 1;
                end
                
                channel = Options.Preprocessing(siteIndex).MaskChannel;
                
                
                for cropIndex = remainingCrops
                    objectIndex = getObjectIndex(file, site, cropIndex, Options.Objects);
                    for channel = Options.Preprocessing(siteIndex).ChannelsToUse %setdiff(Options.Preprocessing(siteIndex).ChannelsToUse, Options.Preprocessing(siteIndex).MaskChannel)
                        
                        fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
                        stack2 = readFileToStack([ExpSettings.NoBackgroundPath, fileKey, '.tiff']);
                        
                        if isempty(Options.Objects(objectIndex).CropArea)
                            continue
                        else
                            cropArea = [max(Options.Objects(objectIndex).CropArea(1), 1), max(Options.Objects(objectIndex).CropArea(2), 1),...
                                min(Options.Crop(siteIndex).OriginalSize(2) - Options.Objects(objectIndex).CropArea(1), Options.Objects(objectIndex).CropArea(3)),...
                                min(Options.Crop(siteIndex).OriginalSize(1) - Options.Objects(objectIndex).CropArea(2), Options.Objects(objectIndex).CropArea(4))];
                        end              

                        % Make sure the cropping area does not exceed the
                        % image size
                        if cropArea(1) + cropArea(3) > size(stack2, 2) || cropArea(2) + cropArea(4) > size(stack2, 1)
                            cropArea(3) = min(size(stack2,2) - cropArea(1), cropArea(3));
                            cropArea(4) = min(size(stack2,1) - cropArea(2), cropArea(4));
                        end

                        if ~isempty(Options.Objects(objectIndex).CropArea) % Manual or polygon cropping
                            cropMask = zeros(size(stack2, [1,2]));
                            cropMask(cropArea(2):cropArea(2)+cropArea(4),...
                                cropArea(1):cropArea(1)+cropArea(3)) = 1;
                        else
                            continue
                        end                        

                        % Make sure the crop mask is of the same size as
                        % the image (is corrected for jitter)
                        if ~isequal(size(stack2, [1,2]), size(cropMask)) && channel == Options.Preprocessing(siteIndex).ChannelsToUse(1)
                            cropMask = imcrop(cropMask, [Options.Preprocessing(siteIndex).FourierJitter, ...
                                Options.Preprocessing(siteIndex).FourierJitter, size(stack2,2) - 1, size(stack2,1) - 1]);                            
                        end

                        % Mask
                        stack2 = stack2.*uint16(repmat(cropMask, [1,1, size(stack2,3)]));                        
                        
                        % Crop
                        cropArea = floor(cropArea); cropArea(cropArea == 0) = 1; cropArea = uint16(cropArea);
                        if size(stack2, 3) > 1
                            stack3 = imcrop3(stack2, [cropArea(1), cropArea(2), 1,...
                                cropArea(3), cropArea(4), size(stack2,3) - 1]);
                        else
                            stack3 = imcrop(stack2, [cropArea(1), cropArea(2),...
                                cropArea(3), cropArea(4)]);
                        end
                        
%                         if strcmp(Options.Crop(siteIndex).CropMode, 'Automatic')
%                             temp = zeros(Options.Objects(objectIndex).maskSize);
%                             temp(Options.Objects(objectIndex).PixelIdxList) = 1;
%                             temp = imcrop(temp, [cropArea(1), cropArea(2),...
%                             cropArea(3), cropArea(4)]);
%                             stack2 = double(temp).*double(stack2);
%                         end                                                

                        % Save results
                        fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                        uiprogressdlg(uifig, 'Title', 'Cell cropping', 'Message', ['Saving cropped file ', fileKey], 'Indeterminate', 'on');                        
                        Stack2TIFF(stack3, [ExpSettings.CroppedPath, filesep, fileKey, '.tiff']);
                        clear stack3;
                    end % For channel                    
                end % For crop
                remainingCrops = [];
                fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, 0);
                uialert(uifig, {'Cropped file '; fileKey}, 'Success!', 'Icon', 'success');
            end % For site
        end % For file                
    else
        ExpSettings.CroppedPath = ExpSettings.NoBackgroundPath;
    end % If crop
    uiprogressdlg(uifig,'Title',"Chop chop! Let's crop!", 'Message', 'Saving ExpSettings and Options data...', 'Indeterminate', 'on');
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-ExpSettings.mat'], 'ExpSettings');
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-Options.mat'], 'Options', '-v7.3');    
    close(uifig);
    warning on
%
% for file = ExpSettings.Selected.files_to_analyze
% for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
% siteIndex = getSiteIndex(file, site, metadata);
% if metadata(file).n_wells > 1
% well = find(site <= cumsum(metadata(file).wells(:,3)), 1);
% else
% well =1;
% end
% Options.Crop(siteIndex).File = Options.Preprocessing(siteIndex).File;
% Options.Crop(siteIndex).Site = Options.Preprocessing(siteIndex).Site;
% remainingCrops = 1:size(Options.Crop(siteIndex).Cropped, 2);
% for cropIndex = remainingCrops
%     cropArea = Options.Crop(siteIndex).Cropped(cropIndex).CropArea;
%     channel = Options.Preprocessing(siteIndex).MaskChannel;
%     fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);          
%     stack1 = readFileToStack([ExpSettings.NoBackgroundPath, fileKey, '.tif']);
% 
%     stack1 = imcrop3(stack1, [cropArea(1), cropArea(2), 1,...
%     cropArea(3), cropArea(4), size(stack1,3)-1]);
% 
%     fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
%     Stack2TIFF(stack1, [ExpSettings.CroppedPath, filesep, fileKey, '.tif']);
% 
%                 
%     for channel = setdiff(Options.Preprocessing(siteIndex).ChannelsToUse, Options.Preprocessing(siteIndex).MaskChannel)
%         fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well);
%         stack2 = readFileToStack([ExpSettings.NoBackgroundPath, fileKey, '.tif']);
%         stack2 = imcrop3(stack2, [cropArea(1), cropArea(2), 1,...
%         cropArea(3), cropArea(4), size(stack2,3)-1]);
%         fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
%         Stack2TIFF(stack2, [ExpSettings.CroppedPath, filesep, fileKey, '.tif']);
%     end
%     clear stack1 stack2;
%     fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
%     disp(fileKey);
% end
% 
% end
% end

end

%% Segmentation (Masking)
for s = 1 % Fake for loop to enable section collapse
    if isfield(Options.Preprocessing, 'MaskChannel')
        if exist(ExpSettings.CroppedPath, 'dir') && isfield(Options, 'Crop')
            sourcePath = ExpSettings.CroppedPath;
        else
            sourcePath = ExpSettings.NoBackgroundPath;
        end
        
        selection = 'None';
        
        uifig = uifigure('Name', 'Segmentation time, come on!'); uifig.WindowStyle = 'alwaysontop';
        supervise = uiconfirm(uifig, 'Do you want to supervise the progress? ','Select mode',...
            'Options', {'Yes', 'No, run them all, I will check them later', 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 3);

        for file = ExpSettings.Selected.files_to_analyze
            if ~exist('remainingSites', 'var') || isempty(remainingSites)
                remainingSites = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file));
            end
%             fg = figure; ax1 = gca; hold(ax1, 'on'); fg.Visible = 0;
                        
            for site = remainingSites
                well = getWellIndex(file, site, metadata);
                siteIndex = getSiteIndex(file, site, metadata);
                
                if isfield(Options, 'Crop')
                    if ~exist('remainingCrops', 'var') || isempty(remainingCrops)
                        remainingCrops = 1:Options.Crop(siteIndex).NumCrops;
                    end                    
                else
                    remainingCrops = 0;
                end

                for cropIndex = remainingCrops
                    if isfield(Options, 'Crop')
                        objectIndex = getObjectIndex(file, site, cropIndex, Options.Objects);
                        if isempty(Options.Objects(objectIndex).CropArea)                            
                            continue                            
                        end
                        
                        if ~Options.Objects(objectIndex).Retained
                            continue
                        end
                    end
                    
                    fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention,well, cropIndex);
                    uiprogressdlg(uifig, 'Title', ['Segmenting file ', fileKey], 'Message', 'Loading data...', 'Indeterminate', 'on');

                    fileKey = getFileKey(file, site, Options.Preprocessing(siteIndex).MaskChannel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                    stack_path = [sourcePath, fileKey, '.tiff'];
                    stack = readTIFFstack(stack_path);
                    stack = double(stack);

                    if strcmpi(supervise, 'Yes')
                        approved = 0;
                        fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention,well, cropIndex);
                        playMovie = uiconfirm(uifig, 'Do you want to watch all the frames masked? ',...
                            ['Play movie ', fileKey],...
                            'Options', {'Yes', 'No, just the last one', 'Cancel'}, 'DefaultOption', 2, 'CancelOption', 3);
                        uifig.WindowStyle = 'alwaysontop';
                        
                        uiprogressdlg(uifig, 'Title', ['Segmenting file ', fileKey], 'Message', 'Masking...', 'Indeterminate', 'on');
                        [~, mask, TempOpt] = getCellMask(stack, Options.Preprocessing, siteIndex, cropIndex);

                        fig = figure; 
                        fig.WindowState = 'maximized';
                        uifig.WindowState = 'minimized';
                        uifig.WindowState = 'normal';

                        ax = subplot(5,1,1:4); hold(ax, 'on'); 
                        title(ax, 'Segmentation preview');
                        subtitle(ax, ['Sharpen Method = ', Options.Preprocessing(siteIndex).SharpenMethod, ...
                            '. Mask Method = ', Options.Preprocessing(siteIndex).MaskMethod]);

                        ax2 = subplot(5,1,5); hold(ax2, 'on');
                        xlim(ax2, [0, size(stack,3)]);
                        title(ax2, 'Segmentation thresholds');
                        xlabel(ax2, 'Frame');
                        ylabel(ax2, 'Intensity');

                        if strcmpi(playMovie, 'Yes')
                            frameRange = 1:size(stack,3);
                        else
                            frameRange = size(stack,3);
                        end

                        plot(ax2, TempOpt(siteIndex).Threshold(:, max(cropIndex, 1)));
                        xl = xline(ax2, 1, 'LineStyle', '--', 'Color', 'r');
                        for frame = frameRange
                            image = stack(:,:,frame);
                            if sum(mask, "all") ~= 0
                                outlined_img = DrawMaskOutline(image, mask(:,:,frame));
                                imagesc(ax, outlined_img); axis(ax, 'image'); axis(ax, 'ij'); 
                                xl.Value = frame;
                                ax.Title.String = ['Segmentation preview. Frame # ', num2str(frame)];                                                                
                                drawnow;
                            else
                                disp(['Warning! Null mask at frame # ', num2str(frame)]);
                            end            
                        end

                        selection = uiconfirm(uifig, 'Accept segmentation? ', 'Automatic segmentation',...
                            'Options', {'Accept', 'Improve', 'Cancel'}, 'DefaultOption', 1, 'CancelOption', 3);
                        uifig.WindowState = 'minimized';
                        uifig.WindowState = 'normal';

                        switch selection
                            case 'Accept'
                                approved = 1;
                                if (length(remainingCrops) == 1 || cropIndex == remainingCrops(end))
                                    remainingSites = ExpSettings.Selected.sites(file, find(ExpSettings.Selected.sites(file,...
                                        1:ExpSettings.Selected.n_sites_per_file(file)) == site)+1:...
                                        ExpSettings.Selected.n_sites_per_file(file));
                                    remainingCrops = [];
                                elseif length(remainingCrops) > 1 && cropIndex ~= remainingCrops(end)
                                    remainingSites = ExpSettings.Selected.sites(file, find(ExpSettings.Selected.sites(file,...
                                        1:ExpSettings.Selected.n_sites_per_file(file)) == site):...
                                        ExpSettings.Selected.n_sites_per_file(file));
                                    remainingCrops = setdiff(remainingCrops, cropIndex);
                                end
%                                 close(fig);
                            case 'Improve'
                                approved = 0;
                                remainingSites = ExpSettings.Selected.sites(file, find(ExpSettings.Selected.sites(file,...
                                    1:ExpSettings.Selected.n_sites_per_file(file)) == site):...
                                    ExpSettings.Selected.n_sites_per_file(file));
                                Options.Preprocessing(siteIndex).Threshold = TempOpt(siteIndex).Threshold;
%                                 close(fig);
                                open('segmentationOptimization.mlx');
                                break
                                % Let the user optimize the segmentation, save the parameters and stop.
                                % The code has to be run again for the selected sites,
                                % this time the optimized parameters will be directly retrieved
                                % from the Options structure.
                            case 'Cancel'
                                approved = 0;
                                close(fig);
                                continue % Skip this site and continue with the rest
                        end
                    elseif ~strcmpi(supervise, 'Cancel')
                        approved = 1;
                        selection = 'Accept';
                    end

                    if approved
                        fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention,well, cropIndex);
                        h = uiprogressdlg(uifig, 'Title', ['Segmenting file ', fileKey], 'Message', 'Masking...', 'Indeterminate', 'on'); 
                        if strcmpi(supervise, 'Yes')
                            SegOpt = TempOpt;
                        else
                            [~, mask, SegOpt] = getCellMask(stack, Options.Preprocessing, siteIndex, cropIndex);
                        end

                        h.Message = 'Saving mask...';
                        fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                        Options.Preprocessing = SegOpt;            
                        save([ExpSettings.MasksPath, fileKey, '-Mask.mat'], 'mask', '-v7.3');

                        h.Message = 'Saving outlined mask TIFF...';
                        outlined_stack = cell(1,size(mask,3));
                        for frame = 1:size(mask,3)
                            if frame == 1; WriteMode = "overwrite"; else; WriteMode = "append"; end
                            imwrite(uint16(mask(:,:,frame)), [ExpSettings.MasksPath, fileKey, '-Mask.tif'],...
                                "WriteMode", WriteMode, "Compression", "none");

                            outlined_stack{frame} = DrawMaskOutline(stack(:,:,frame), mask(:,:,frame), [0,1,0]);
                            imwrite(outlined_stack{frame}, ...
                                [ExpSettings.MasksPath, fileKey,'-Outlined-Stack.tif'], ...
                                'WriteMode', WriteMode, 'Compression', 'none');                            
                        end                            
                        clear outlined_stack;

                        % h.Message = 'Saving masked stack...';
                        % for channel = Options.Preprocessing(siteIndex).ChannelsToUse
                        %     % Mask stack
                        %     fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);            
                        %     stack_path = [sourcePath, fileKey, '.tiff'];
                        %     stack = readTIFFstack(stack_path);
                        %     masked_stack = double(stack).*mask; masked_stack(masked_stack == 0) = NaN;          
                        %     save([ExpSettings.MasksPath, fileKey, '-Masked-Stack.mat'], 'masked_stack', '-v7.3');
                        % end                    

                        clear mask; %masked_stack;                
                    end

%                         plot(ax1, Options.Preprocessing(siteIndex).Threshold);
%                         hold(ax1, 'on');

                    Options.Preprocessing(siteIndex).DateModified = datestr(now);
                    clear stack;
                    fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                end
                
                switch selection
                    case 'Accept'
                        if length(remainingCrops) == 1
                            remainingSites = ExpSettings.Selected.sites(file, find(ExpSettings.Selected.sites(file,...
                                1:ExpSettings.Selected.n_sites_per_file(file)) == site)+1:...
                                ExpSettings.Selected.n_sites_per_file(file));
                        end
                        remainingCrops = [];
                        fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention);
                        uialert(uifig, 'Success!', ['File ', fileKey], 'Icon', 'success');
                        pause(1);
                    case 'Improve'
                        remainingSites = ExpSettings.Selected.sites(file, find(ExpSettings.Selected.sites(file,...
                            1:ExpSettings.Selected.n_sites_per_file(file)) == site):...
                            ExpSettings.Selected.n_sites_per_file(file));
                        open('segmentationOptimization.mlx');
                        break
                        % Let the user optimize the segmentation, save the parameters and stop.
                        % The code has to be run again for the selected sites,
                        % this time the optimized parameters will be directly retrieved
                        % from the Options structure.
                    case 'Cancel'                        
                        continue % Skip this site and continue with the rest
                end

            end

%             fg.Visible = 1;
%             savefig(fg, [ExpSettings.MasksPath, ExpSettings.ExpID, '-', num2str(file), '-MaskThresholdsPlot.fig']);
            
        end        
        
        close(uifig);
        save([ExpSettings.CloudPath, ExpSettings.ExpID, '-ExpSettings.mat'], 'ExpSettings');
        save([ExpSettings.CloudPath, ExpSettings.ExpID, '-Options.mat'], 'Options', '-v7.3');
    end
end

%% Possible Next Steps:
% Generate FRET maps:
open FRET.m

% Extract Objects and their properties:
open ObjectProperties.m

%% Sharpening (for STICS on myosin)
for h = 1 % Fake for loop to enable section collapse
    
    if exist(ExpSettings.CroppedPath, 'dir')
        sourcePath = ExpSettings.CroppedPath;
    else
        sourcePath = ExpSettings.NoBackgroundPath;
    end
    ExpSettings.SharpPath = [sourcePath, 'Sharpenned', filesep];
    if ~exist(ExpSettings.SharpPath, 'dir'); mkdir(ExpSettings.SharpPath); end
        
    for file = ExpSettings.Selected.files_to_analyze
        for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
            siteIndex = getSiteIndex(file, site, metadata);
            Options.Preprocessing(siteIndex).SharpGaussianSize = 3;
            
            channel = 3;
            
            if strcmpi(sourcePath, ExpSettings.CroppedPath)
                remainingCrops = 1:size(Options.Crop(siteIndex).Cropped, 2);
            else
                remainingCrops = 1;
            end

            for cropIndex = remainingCrops                    
                fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                disp(fileKey);
                stack_path = [sourcePath, fileKey, '.tif'];
                stack = readFileToStack(stack_path);
                
                GaussianBlur = imgaussfilt3(stack, Options.Preprocessing(siteIndex).SharpGaussianSize);
                sharpStack = (stack - GaussianBlur);
                
                Stack2TIFF(sharpStack, [ExpSettings.SharpPath, fileKey, '-SH.tiff']);
            end
            
        end
    end
    
    save([ExpSettings.CloudPath, ExpSettings.ExpID, '-ExpSettings.mat'], 'ExpSettings');
    save([ExpSettings.CloudPath, ExpSettings.ExpID, '-Options.mat'], 'Options', '-v7.3');    
    
end

%% Heterogeneity Maps
for h = 1 % Fake for loop to enable section collapse  
    if exist(ExpSettings.CroppedPath, 'dir')
        sourcePath = ExpSettings.CroppedPath;
    else
        sourcePath = ExpSettings.NoBackgroundPath;
    end
    ExpSettings.HetMapPath = [sourcePath, 'Heterogeneity Maps', filesep];
    if ~exist(ExpSettings.HetMapPath, 'dir'); mkdir(ExpSettings.HetMapPath); end
    
    for file = ExpSettings.Selected.files_to_analyze
        for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
            siteIndex = getSiteIndex(file, site, metadata);
            
            if metadata(file).n_wells > 1
                well = find(site <= cumsum(metadata(file).wells(:,3)), 1);
            else
                well =1;
            end
            
            Options.Preprocessing(siteIndex).HetMapROIsize = 8;
            
            channel = 1;
            
            if strcmpi(sourcePath, ExpSettings.CroppedPath)
                nCrops = 1:Options.Crop(siteIndex).NumCrops;
            else
                nCrops = 1;
            end

            for cropIndex = nCrops
                
                fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                load([ExpSettings.MasksPath, fileKey, '-Mask.mat'], 'mask');
            
                fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                disp(fileKey);
                stack_path = [sourcePath, fileKey, '.tiff'];
                stack = readFileToStack(stack_path);
                
                heterogeneityMap = getHeterogeneityMap(stack, Options.Preprocessing(siteIndex).HetMapROIsize);
%                 deltaHetMap = diff(heterogeneityMap, 1, 3);
%                 heterogeneityMap = heterogeneityMap + 10000;
                heterogeneityMap = double2uint(heterogeneityMap, 8, [-0.5, 2]); % For HM using sqrt(std);
%                 heterogeneityMap = double2uint(heterogeneityMap, 8, [0, 255]); % For HM using STD
                
%                 heterogeneityMap(heterogeneityMap == 2^16) = 0;
                heterogeneityMap(~mask) = 0;
                
                Stack2TIFF(heterogeneityMap, [ExpSettings.HetMapPath, fileKey, '-HMD-',...
                    num2str(Options.Preprocessing(siteIndex).HetMapROIsize), '.tif']);                                
                                
            end                        
        end
    end
end

%% Create STICS Options container
for o = 1 % Fake for loop to enable section collapse
    
    if exist(ExpSettings.SharpPath, 'dir')
        sourcePath = ExpSettings.SharpPath;
    else
        sourcePath = ExpSettings.CroppedPath;
    end
    ExpSettings.STICS_Path = [ExpSettings.RootPath, 'STICS', filesep];
    if ~exist(ExpSettings.SharpPath, 'dir'); mkdir(ExpSettings.SharpPath); end
    
    index = 0;
    for file = ExpSettings.Selected.files_to_analyze
        for site = 1:metadata(file).n_sites
            siteIndex = getSiteIndex(file, site, metadata);

            channel = 3;
            
            if exist(ExpSettings.CroppedPath, 'dir')
                remainingCrops = 1:size(Options.Crop(siteIndex).Cropped, 2);
            else
                remainingCrops = 1;
            end

            for cropIndex = remainingCrops                    
                index = index +1;
                fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention, well, cropIndex);
                
                
                Options.STICS(index).File = file;
                Options.STICS(index).Site = site;
                Options.STICS(index).Crop = cropIndex;
                Options.STICS(index).FileKey = fileKey;
                Options.STICS(index).Settings = struct();
                
                disp(fileKey);
            end
            
        end
    end
    
    save([ExpSettings.CloudPath, ExpSettings.ExpID, '-ExpSettings.mat'], 'ExpSettings');
    save([ExpSettings.CloudPath, ExpSettings.ExpID, '-Options.mat'], 'Options', '-v7.3');    
    
end

%% Photobleach correction for TIFF Files (not FRET)
for p = 1 % Fake for loop to enable section collapse
%     correctBleaching(ExpSettings, Options.Preprocessing, metadata);
    if ~isfield(ExpSettings, 'BleachCorrected'), ExpSettings.BleachCorrected = [sourcePath, 'Bleach corrected', filesep]; end
    if ~exist(ExpSettings.BleachCorrected, "dir"); mkdir(ExpSettings.BleachCorrected); end    
    
    for file = ExpSettings.Selected.files_to_analyze        
        combine_sites = 1; % By default
        for channel = Options.Preprocessing(siteIndex).ChannelsToUse
            all_sites_normalized_masked_mean = nan(metadata(file).n_frames, metadata(file).n_sites);            
            for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
                % Select file
                siteIndex = getSiteIndex(file, site, metadata);                
                fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention);            
                disp(fileKey);
                load([ExpSettings.MasksPath, fileKey, '-Masked-Stack.mat'], 'masked_stack');                               
                
                % Estimate photobleaching. 
                % Note: `masked_mean` is the equivalent of `bleach_raw`.
                [masked_mean, fitpara, corr] = estimatePhotobleaching(masked_stack, 'mean',...
                    Options.Preprocessing(siteIndex).BleachingModel, 1);
                save([ExpSettings.MasksPath, fileKey, '-Masked-Mean.mat'], 'masked_mean');
                
                normalized_masked_mean = masked_mean./median(masked_mean, 1, 'omitnan');                 
                save([ExpSettings.MasksPath, fileKey, '-Normalized-Masked-Mean.mat'], 'normalized_masked_mean');
                all_sites_normalized_masked_mean(:, site) = normalized_masked_mean; % Equivalent to `bleach_raw_all`.                
                

                combine_sites = logical(Options.Preprocessing(siteIndex).BleachCorrectPerFile * combine_sites);
                
                if ~Options.Preprocessing(siteIndex).BleachCorrectPerFile || ~combine_sites
                    if ~isfield(Options.Preprocessing(siteIndex), 'BleachParameters') || size(Options.Preprocessing(siteIndex).BleachParameters, 1) == metadata(file).n_channels
                        Options.Preprocessing(siteIndex).BleachParameters = {};
                        c_index = 1;
                    else
                        c_index = size(Options.Preprocessing(siteIndex).BleachParameters, 1) + 1;
                    end
                    Options.Preprocessing(siteIndex).BleachParameters{c_index, 1} = fitpara;

                    % Load stack
                    stack_path = [sourcePath, fileKey, '.tif'];
                    stack = readTIFFstack(stack_path);
                    
                    % Load mask
                    fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention);
                    load([ExpSettings.MasksPath, fileKey, '-Mask.mat'], 'mask');
                    
                    % Correct single site
                    corr_norm = corr./median(corr);
                    bc_stack = zeros(size(stack));
                    mean_bc_stack = zeros(size(stack,3), 1);                    
                    fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention);
                    for frame = 1:metadata(file).n_frames
                        bc_stack(:,:,frame) = double(stack(:,:,frame)./(corr_norm(frame, 1)));
                        masked_bc_stack = (bc_stack(:,:,frame).*mask(:,:,frame));
                        masked_bc_stack(masked_bc_stack == 0) = NaN;
                        mean_bc_stack(frame) = mean(masked_bc_stack, "all", "omitnan");
                        if frame == 1; WriteMode = "overwrite"; else; WriteMode = "append"; end
                        imwrite(uint16(bc_stack(:,:,frame)), [ExpSettings.BleachCorrected, fileKey, '.tif'],...
                            "WriteMode", WriteMode, "Compression", "none");
                    end                    
                    
                    Options.Preprocessing(siteIndex).DateModified = datestr(now);
                    clear bc_stack stack;
                    
                    % Display results for single site. 
                    if Options.Preprocessing(siteIndex).PlotBleaching
                        [fig, ax] = plotBleachingFit(normalized_masked_mean, corr_norm);
                        hold(ax, 'on'); 
                        ylabel(ax, 'Normalized mean intensity of the masked area');
                        ttl = ax.Title; ttl.String = [ttl.String, ' for site ' fileKey]; ax.Title = ttl;
                        plt = plot(ax, 1:size(mean_bc_stack,1), (mean_bc_stack./median(mean_bc_stack, 1, "omitnan")), 'k');
                        plt.DisplayName = ['Bleach-corrected (', Options.Preprocessing(siteIndex).BleachingModel, ' model)'];
                        drawnow;
                        if Options.Preprocessing(siteIndex).SaveBleachingPlot
                            exportgraphics(ax, [ExpSettings.BleachCorrected, fileKey, '-Bleaching-Profiles.png']);
                        end
                    end                                        
                end

                Options.Preprocessing(siteIndex).DateModified = datestr(now);
            end                                                                                   
            
            if combine_sites
                [median_all_sites_normalized_masked_mean, fitpara, corr] = estimatePhotobleaching(all_sites_normalized_masked_mean, ...
                    'median', Options.Preprocessing(siteIndex).BleachingModel, 0);
                % Note: `median_all_sites_normalized_masked_mean` = `bleach_mean`
                corr_norm = corr./median(corr);
                
                PlotBleaching = 1;
                for site = ExpSettings.Selected.sites(file, 1:ExpSettings.Selected.n_sites_per_file(file))
                    siteIndex = getSiteIndex(file, site, metadata);
                    if ~isfield(Options.Preprocessing(siteIndex), 'BleachParameters') || size(Options.Preprocessing(siteIndex).BleachParameters, 1) == metadata(file).n_channels
                        Options.Preprocessing(siteIndex).BleachParameters = {};
                        c_index = 1;
                    else
                        c_index = size(Options.Preprocessing(siteIndex).BleachParameters, 1) + 1;
                    end
                    Options.Preprocessing(siteIndex).BleachParameters{c_index, 1} = fitpara;                    

                    fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention);
                    
                    % Load stack
                    stack_path = [sourcePath, fileKey, '.tif'];
                    stack = readTIFFstack(stack_path);
                    
                    % Load mask
                    fileKey = getFileKey(file, site, 0, metadata, ExpSettings.Selected.Convention);
                    load([ExpSettings.MasksPath, fileKey, '-Mask.mat'], 'mask');

                    % Apply bleach correction
                    bc_stack = zeros(size(stack));
                    mean_bc_stack = zeros(size(stack,3), 1);                    
                    fileKey = getFileKey(file, site, channel, metadata, ExpSettings.Selected.Convention);
                    for frame = 1:metadata(file).n_frames
                        bc_stack(:,:,frame) = double(stack(:,:,frame)./(corr_norm(frame, 1)));
                        masked_bc_stack = (bc_stack(:,:,frame).*mask(:,:,frame));
                        masked_bc_stack(masked_bc_stack == 0) = NaN;
                        mean_bc_stack(frame) = mean(masked_bc_stack, "all", "omitnan");
                        if frame == 1; WriteMode = "overwrite"; else; WriteMode = "append"; end
                        imwrite(uint16(bc_stack(:,:,frame)), [ExpSettings.BleachCorrected, fileKey, '.tif'],...
                            "WriteMode", WriteMode, "Compression", "none");
                    end
                    
                    PlotBleaching = PlotBleaching * Options.Preprocessing(siteIndex).PlotBleaching;
                    
                    Options.Preprocessing(siteIndex).DateModified = datestr(now);
                    clear bc_stack stack;
                end
                
                % Display result for multiple sites.
                if PlotBleaching
                    [fig, ax] = plotBleachingFit(all_sites_normalized_masked_mean, corr);
                    hold(ax, 'on'); 
                    ylabel(ax, 'Normalized mean intensity of the masked area');
                    ttl = ax.Title; ttl.String = [ttl.String, ' for site ' fileKey]; ax.Title = ttl;
                    plt = plot(ax, 1:size(mean_bc_stack,1), (mean_bc_stack./median(mean_bc_stack, 1, "omitnan")), 'k');
                    plt.DisplayName = ['Bleach-corrected (', Options.Preprocessing(siteIndex).BleachingModel, ' model)'];
                    drawnow;
                    if Options.Preprocessing(siteIndex).SaveBleachingPlot
                        exportgraphics(ax, [ExpSettings.BleachCorrected, fileKey, '-Bleaching-Profiles.png']);
                    end
                end
                
            end
        end        
    end
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-Options.mat'], 'Options');
    clear masked_stack mask stack masked_bc_stack;        
    save([ExpSettings.CloudPath, ExpSettings.ExpID,  '-ExpSettings.mat'], 'ExpSettings');
end

%% Small compatibility bypass if doing STICS:
    if size(mask,3)>1 && isa(mask, 'logical')
        maskFinal = cell(1, size(mask,3));
        for frame = 1:size(mask,3)
            maskFinal{1,frame} = mask(:,:,frame);
        end
    end
%%
    for frame = 1:size(maskFinal,2)
        mask(:,:,frame) = maskFinal{1, frame};
    end
    
%%
% %% %%%% Bleaching correction: Detrmine linear fit parameters for FRET/CFP decay
%             % bleach_1=nanmean(vect(imRatio_raw{1}));
%             for frame=1:length(imRatio_raw)
%                bleach_raw(frame)=nanmean(vect(imRatio_raw{frame}));
%                bleach_raw_mRuby(frame)= nanmean(vect(im_mRub3_raw{frame}));
%             end
% 
%             save([ExpSettings.NoBackgroundPath, fileKey,'_Bleach_raw.mat'],'bleach_raw', 'bleach_raw_mRuby');
%             clear('cellCoors','imRatio_raw', 'bleach_raw', 'bleach_raw_mRuby', 'im_mRub3_raw');
%             for frame = 1:length(maskFinal)
%                 imFRETOutline{frame}=DrawMaskOutline(imFRET_fr{frame},maskFinal{frame});
%             end
%             save([ExpSettings.NoBackgroundPath, fileKey, '_RatioData_raw.mat'],'imFRETOutline','maskFinal', '-append');
%             clear('imFRETOutline', 'maskFinal');        
% 
% %% %% Test parameters using a single frame
% 
% frame = 15;
% 
% % Load raw images
% fileKey = getFileKey(file, site, 1, metadata, Convention);
% imCH1_raw = double(imread([ExpSettings.TIFFStacksPath, fileKey, '.tif'], frame));
% fileKey = getFileKey(file, site, 2, metadata, Convention);
% imCH2_raw = double(imread([ExpSettings.TIFFStacksPath, fileKey, '.tif'], frame));
% 
% % Align raw images
% imstack(:,:,1) = imCH1_raw; imstack(:,:,2) = imCH2_raw; 
% imaligned = dualviewAlignFromFittedSurface(imstack,pX,pY,1);
% imaligned = makeSquare(imaligned);
% imCH1 = imaligned(:,:,1);
% imCH2 = imaligned(:,:,2);
% 
% % Generate background mask, background-subtract raw images
% bgmask = getBGMask(imCH1+imCH2);
% % bgmask = getCellMaskCyto_3(imCH2+imCH1,0.65,7,0,0);
% imCH1bg = subBG(imCH1, bgmask, CH1_AVG_BG_Aligned);
% imCH2bg = subBG(imCH2, bgmask, CH2_AVG_BG_Aligned);
% 
% 
% % Display background mask
% figure;
% subplot(1,2,1); imagesc(bgmask);
% subplot(1,2,2); imagesc(imCH1bg+imCH2bg);
% 
% % Segment cells (generate foreground/background mask)
% [mask_raw, cellCoorsTemp] = getCellMaskCyto_3(imCH1+imCH2,0.65,7,0,0);
% mask = imfilter(mask_raw,fspecial('disk',2),'symmetric');
% 
% % Detrmine ratio
% imCH1bg(~mask) = nan;
% imCH1_filtered = ndnanfilter(imCH1bg,fspecial('disk',3),'replicate');
% imCH2bg(~mask) = nan;
% imCH2_filtered = ndnanfilter(imCH2bg,fspecial('disk',3),'replicate');
% imRatio=imCH1_filtered./imCH2_filtered;
% imRatio(~mask)=nan; 
% imRatioScaled(~mask)=0;
% 
% colorRange = [round(prctile(imRatio(:),2),1),round(prctile(imRatio(:),98),1)];
% tempRATIO = ratio2RGB(imRatio,colorRange);
% 
% % Display results
% figure; subplot(2,2,1);imshow(mat2gray(imCH1+imCH2),[0 0.2]);
% subplot(2,2,2); imshow(mask);
% subplot(2,2,3); imshow(tempRATIO); 
% subplot(2,2,4); imshow(DrawMaskOutline(imCH1+imCH2,mask,[0 1 0]));    
% 
% figure; imagesc(imRatio); axis image; colorbar;
% 
    