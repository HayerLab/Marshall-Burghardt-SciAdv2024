%% SETUP
%% Generate ExpSettings
ExpSettings = generateExpSettings();
save([ExpSettings.CloudPath, ExpSettings.ExpID, '-ExpSettings']);
%% If working with TIFFs (no ND2 files available) 
% Generate the metadata manually

metadata = struct();
% metadata.wells(:,1:2) = metadata.wells(:, 1:2) + 1;
%% Move all the files to the same folder:
filelist = dir(fullfile(ExpSettings.DataPath, '**\*.*'));  %get list of files and folders in any subfolder
filelist = filelist(~[filelist.isdir]);  %remove folders from list

for n = 1:length(filelist)
    expression = '\w*tif';
    matchStr = regexp(filelist(n).name, expression, 'match');
     if strcmp(matchStr, 'tif')
        tobemoved = [filelist(n).folder, filesep, filelist(n).name];
        copyfile(tobemoved, ExpSettings.TIFFStacksPath);
     end
end

%% Rename well names
metadata.wellNames = {};
rows = ['A'; 'B'; 'C'; 'D'; 'E';'F'; 'G'; 'H'];
cols = [1; 2; 3; 4; 5; 6; 7; 8; 9;10; 11; 12];
rowNames = {'','CNTRL', 'CPD31', 'ROCKi', 'NCDZL', 'NOC31', 'NOCY2'}; % Row names
colNames = {'', '', ''}; % Column headers
i = 1;
for row = 2:7
    for col = 2:3
        if ~isempty(colNames{col})
            metadata.wellNames{i,1} = [[rows(row), num2str(cols(col))], '-', rowNames{row}, '-', colNames{col}];
        else
            metadata.wellNames{i,1} = [[rows(row), num2str(cols(col))], '-', rowNames{row}];
        end
        i = i + 1;
    end
end

save([ExpSettings.MetadataPath, ExpSettings.ExpID,'-Metadata.mat'], 'metadata');

%% Change convention
file = 1;%ExpSettings.Selected.files_to_analyze;

filelist2 = dir(fullfile([ExpSettings.DataPath], '**\*.*'));  % get list of files and folders in any subfolder
filelist2 = filelist2(~[filelist2.isdir]); % remove folders from list

for i = 1:length(filelist2)
    % Skip any other filetypes other than TIF or TIFFs
    if ~contains(filelist2(i).name, '.tif')
        continue
    end

    % Extract row, column and site from filename:
    filename = filelist2(i).name;
    expression = '\_'; splitStr = regexp(filename,expression,'split')';
    row = splitStr(1); row = str2num([row{:}]);
    col = splitStr(2); col = str2num([col{:}]);
    siteIndex = splitStr(3); site = str2num(siteIndex{1}); 
    site = getSiteFromRowCol(file, row, col, site, metadata);    
    wellIndex = getWellIndex(file, site, metadata);
    siteIndex = str2num([siteIndex{:}]);
    if siteIndex >= 10
        siteIndex = int2str(siteIndex);
    else
        siteIndex = strcat('0',int2str(siteIndex));
    end

    % Determine channels
    if length(filelist2) == metadata(file).seq_count % One TIFF per channel
        % The channel name is contained in the filename
        channel_exp = '\ ';
        channelName = regexp(splitStr(4), channel_exp, 'split');
        vec = channelName{1,1}';
        x = vec{1}; 
        switch x
            case '640'
                channelName = 'CY5-T';
                channel = 1; %find(metadata(file).channels);
            case '508'
                channelName = 'YFP';
            case '395'
                channelName = 'HOECHST';
                channel = 2;
        end

    % Change filename
    new_filename = getFileKey(file, site, channel, metadata, "R", wellIndex);
    oldfilepath = [filelist2(i).folder, filesep, filename];
    renameFile(oldfilepath, new_filename);
    disp([filename, ' -> ', new_filename]);
        
    elseif mod(length(filelist2), metadata(file).n_channels) == 0
        % Each TIFF contains all channels as layers
        % Assuming layers are in the same order as metadata.channels
        stack = readFileToStack([filelist2(i).folder, filesep, filelist2(i).name]);
        if size(stack, 3) == metadata(file).n_channels
            for channel = 1:metadata(file).n_channels
                new_filename = getFileKey(file, site, channel, metadata, "R", wellIndex);
                % if channel ~= metadata(file).n_channels % First n-1 channels
                    imwrite(stack(:,:,channel), [ExpSettings.TIFFStacksPath, new_filename, '.tiff'],...
                        'Compression', 'none');
                % else % Last channel: overwrite the file and then rename it
                %     oldfilepath = [filelist2(i).folder, filesep, filelist2(i).name];
                %     imwrite(stack(:,:,channel), oldfilepath,'Compression', 'none');
                %     renameFile(oldfilepath, new_filename);                      
                % end
            end
            disp([filelist2(i).name, ' -> ', new_filename]);
        end
    end        
end

%% Rename all files
% folderPath = 'E:\230416\TIFF Stacks\230416-EP-01\';
% prefix = '230416-EP-01-';
% suffix = '';
% 
% index = 1;
% for row = 2:5
%     for col = 7:8
%         regExpCellArray{index, 1} = [num2str(row), '_', num2str(col)];
%         index = index + 1;
%     end
% end
% % regExpCellArray = {'2_7'; '2_8'; '3_7';'};
% replaceWithArray = metadata(3).wellNames;
% for i = 1:length(regExpCellArray)    
%     renameAllFiles(folderPath, regExpCellArray{i},...
%         replaceWithArray{i}, prefix, suffix);
% end
% 
% for site = 1:metadata(3).n_sites
%     renameAllFiles(folderPath, ['_', num2str(site), '_'],...
%         ['-', num2str(site), '-'], '', suffix);
% end
% 
% renameAllFiles(folderPath, '508 nm_1', 'YFP', '', '');

%% Metadata
% Load or select metadata file
if ~exist('metadata', 'var')
    m = getFilenames(ExpSettings.MetadataPath, 'etadata'); 
    % intentionally left without the m so it finds it regardless of the
    % uppercase or lowercase m.
    if size(m,1) == 1
        try
            load([metadata_path, m{1}], "-mat", 'metadata');
        catch
            disp("Couldn't load metadata file automatically, please select it manually.");
            [metadata_filename, metadata_path] = uigetfile('*.mat');
            load(fullfile(metadata_path, metadata_filename), "-mat", 'metadata');
            ExpSettings.MetadataPath = metadata_path;
        end
    elseif size(m,1) == 0

        disp('No metadata file found, extracting metadata automatically...');

        % For ND2 data:
        % Request from user
    %     fileNames = getFilenames(ExpSettings.DataPath, ".nd2");
    %     n_files = size(fileNames,1);
    %     frame_interval = zeros(n_files, 1);
    %     for f = 1:n_files
    %         fileNames{f}(1:end-4)
    %         frame_interval(f,1) = input(['Frame rate (s) in ', fileNames{f}(1:end-4), '? ']);        
    %     end

        % Extract the rest
        addpath(genpath('D:\GitHub\QLS-PhD\External\bfmatlab'));
        metadata = getMetadata(ExpSettings.DataPath);
    end
end

save([ExpSettings.MetadataPath, ExpSettings.ExpID,'-Metadata.mat'], 'metadata');

%% Background Metadata
% For background:
try
    bgm = getFilenames(ExpSettings.BGPath, 'etadata'); 
    if size(bgm,1) == 1
        try
            load([ExpSettings.BGPath, bgm{1}], "-mat", 'BGmetadata');
            disp(['Loaded ', ExpSettings.BGPath, bgm{1}]);
        catch
            disp("Couldn't load background metadata file automatically, please select it manually.");
            [bgmetadata_filename, bgmetadata_path] = uigetfile('*.mat');
            load(fullfile(bgmetadata_path, bgmetadata_filename), "-mat", 'BGmetadata');    
        end
    elseif size(m,1) == 0

        disp('No metadata file found, extracting metadata automatically...');

        % For ND2 data:
        % Request from user
    %     fileNames = getFilenames(ExpSettings.DataPath, ".nd2");
    %     n_files = size(fileNames,1);
    %     frame_interval = zeros(n_files, 1);
    %     for f = 1:n_files
    %         fileNames{f}(1:end-4)
    %         frame_interval(f,1) = input(['Frame rate (s) in ', fileNames{f}(1:end-4), '? ']);        
    %     end

        % Extract the rest
        BGmetadata = getMetadata(ExpSettings.BGPath);
    else
        disp("Too many files found, please the background metadata file manually.");
        [bgmetadata_filename, bgmetadata_path] = uigetfile('*.mat');
        load(fullfile(bgmetadata_path, bgmetadata_filename), "-mat", 'BGmetadata');    
    end    
catch
    fileNames = getFilenames(ExpSettings.BGPath, ".nd2");
    n_files = size(fileNames,1);

    addpath(genpath('D:\GitHub\QLS-PhD\External\bfmatlab'));
    BGmetadata = getMetadata(ExpSettings.BGPath);    
end
save([ExpSettings.BGPath, ExpSettings.ExpID,'-BGMetadata.mat'], 'BGmetadata');

%% Select files
if ~isfield(ExpSettings.Selected, 'sites')
    ExpSettings = selectWhich(metadata, ExpSettings, 'custom');
end
save([ExpSettings.CloudPath, ExpSettings.ExpID, '-ExpSettings']);

%% Next Step: Preprocessing
open Preprocessing.m;
