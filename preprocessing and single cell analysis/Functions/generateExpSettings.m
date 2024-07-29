function ExpSettings = generateExpSettings(rootPath)
if nargin == 1
    rootProvided = 1;
    root = rootPath;
    ExpSettings.RootPath = root;
    addpath(genpath(root));
else
    rootProvided = 0;
end

%% Set Directories
if ~rootProvided
    % Root
    [root] = strcat(uigetdir("*", ...
        "Please select the master folder (root) that contains all the data..." + ...
        " and will contain the heavy results."), ...
        filesep);
    addpath(genpath(root));
    ExpSettings.RootPath = root;
end

% Data (ND2 or TIFF)
[data_path] = strcat(uigetdir(root,...
    "Please select the folder that contains the file(s) to analyse."), ...
    filesep);
ExpSettings.DataPath = data_path;

% Background
[bg_path] = strcat(uigetdir(root, ...
    "Please select the folder that contains the background file(s)."), ...
    filesep);
ExpSettings.BGPath = bg_path;

% Camera alignment
ExpSettings.AlignmentPath = [ExpSettings.RootPath, 'Alignment', filesep];
if ~exist('ExpSettings.AlignmentPath', 'dir')
    mkdir(ExpSettings.AlignmentPath)    
end

% Cloud
try
    OneDrive = strcat(uigetdir('D:\OneDrive - McGill University\Data\Hayer', ...
    "Please select the cloud folder that will contain the results and lightweight files (ex. OneDrive)."), ...
    filesep);
catch
    OneDrive = ExpSettings.RootPath;
end
addpath(genpath(OneDrive));
ExpSettings.CloudPath = OneDrive;


%% Set Convention
convention = input('Convention? Type "R" or "S"'); %type with quotations
ExpSettings.Selected.Convention = convention;


%% Define paths
if strcmpi(convention, "R")
    ExpName = split(root,'\');
    ExpName = ExpName{size(ExpName, 1)-1,1};
    Exp_ID = ExpName(1:6);
    ExpSettings.ExpID = Exp_ID;
end

ExpSettings.TIFFStacksPath = [ExpSettings.RootPath, 'TIFF Stacks', filesep];
ExpSettings.AlignmentPath = [ExpSettings.RootPath, 'Alignment', filesep];
ExpSettings.BGPath = [ExpSettings.RootPath, 'Background', filesep];
ExpSettings.JitterPath = [ExpSettings.RootPath, 'Jitter', filesep];

if strcmpi(ExpSettings.CloudPath, ExpSettings.RootPath)
    masks_path = [ExpSettings.CloudPath, 'Masks', filesep];
    if ~exist(masks_path, "dir")
        mkdir(masks_path);
    end
    fret_path = [ExpSettings.CloudPath, 'FRET', filesep];
    if ~exist(fret_path, "dir")
        mkdir(fret_path);
    end
        
else
    masks_path = [ExpSettings.RootPath, 'Masks', filesep];
    if ~exist(masks_path, "dir")
        mkdir(masks_path);
    end
    fret_path = [ExpSettings.RootPath, 'FRET', filesep];
    if ~exist(fret_path, "dir")
        mkdir(fret_path);
    end    
end
ExpSettings.MasksPath = masks_path;
ExpSettings.FRET = fret_path;

% paths = fieldnames(ExpSettings);
% paths = paths(boolRegExp(paths, 'Path'));
% for f = 1:size(paths,1)
%     if ~exist(['ExpSettings.', paths{f,1}], 'dir')
%         mkdir(['ExpSettings.', paths{f,1}]);
%     end
% end

%% Metadata
metadata_path = [ExpSettings.CloudPath, 'Metadata', filesep];
if ~exist('metadata_path', 'dir')
    mkdir(metadata_path);
end

addpath(genpath(metadata_path));
ExpSettings.MetadataPath = metadata_path;


end