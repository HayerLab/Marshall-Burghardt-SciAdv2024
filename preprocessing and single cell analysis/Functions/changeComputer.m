function ExpSettings = changeComputer(ExpSettings, mode)
    if nargin < 2
        mode = 'default';
    end    
    originalRootFolder = ExpSettings.("RootPath");
    if ~isfolder(ExpSettings.RootPath) || strcmpi(mode, 'force')
        uifig = uifigure;
        uialert(uifig, 'Please select the corresponding root folder in this computer.',...
            'Directory not found.');
        pause(3);
        newRootFolder = uigetdir(path, 'Please select the corresponding root folder in this computer.');
        newRootFolder = [newRootFolder, filesep];
    else
        newRootFolder = ExpSettings.RootPath;
    end

    originalCloudFolder = ExpSettings.("CloudPath");
    if ~isfolder(ExpSettings.CloudPath) || strcmpi(mode, 'force')
        if ~ishandle(uifig); uifig = uifigure; end
        uialert(uifig, 'Please select the corresponding cloud folder in this computer.',...
            'Directory not found.');
        pause(3);
        newCloudFolder = uigetdir(path, 'Please select the corresponding cloud folder in this computer.');    
        newCloudFolder = [newCloudFolder, filesep];        
    else
        newCloudFolder = ExpSettings.CloudPath;
    end

    fNames = fieldnames(ExpSettings);
        for f = 1:length(fNames)
            value = getfield(ExpSettings, fNames{f});
            if ischar(value) && contains(value, originalRootFolder)
                newValue = replace(value, originalRootFolder, newRootFolder);
                ExpSettings = setfield(ExpSettings,fNames{f},newValue);
            elseif ischar(value) && contains(value, originalCloudFolder)
                newValue = replace(value, originalCloudFolder, newCloudFolder);
                ExpSettings = setfield(ExpSettings,fNames{f},newValue);
            else
                continue
            end
        end

    fNames = fieldnames(ExpSettings);
        for f = 1:length(fNames)
            value = getfield(ExpSettings, fNames{f});
            if ischar(value) && ~isfolder(value)
                mkdir(value)
            end
        end
end
