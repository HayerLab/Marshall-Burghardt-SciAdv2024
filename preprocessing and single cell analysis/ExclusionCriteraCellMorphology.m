%% Exclusion Criteria

% Note that only Options.Objects.Retained is modified. Objects is untouched
% to avoid data loss.

for p = 1:length(Options.Objects)
    Options.Objects(p).Retained = 1;
end

% Find sites with no object found
missingObjectsList = Options.Objects(arrayfun(@(x) isempty(x.CropArea), Options.Objects));
missingIndexes = find(arrayfun(@(x) isempty(x.CropArea), Options.Objects));
for m = 1:length(missingIndexes)
    Options.Objects(missingIndexes(m)).Retained = 0;
end

% Find border objects
borderObjectsList = Options.Objects(arrayfun(@(x) x.BoundaryCell == 1, Options.Objects));
borderIndexes = find(arrayfun(@(x) x.BoundaryCell == 1, Options.Objects));
for b = 1:length(borderIndexes) % Label them as excluded.
    Options.Objects(borderIndexes(b)).Retained = 0;
end

% Find objects detected for cropping, but lost during object extraction
lostObjectsList = Options.Objects(arrayfun(@(x) x.NumObjects2D == 0, Options.Objects));
lostIndexes = find(arrayfun(@(x) x.NumObjects2D == 0, Options.Objects));
for l = 1:length(lostIndexes) % Label them as excluded.
    Options.Objects(lostIndexes(l)).Retained = 0;
end

% Find objects that contain more than one cell touching or bi-nucleate cells
binucleateObjectsList = Options.Objects(arrayfun(@(x) x.NumNuclei ~= 1, Options.Objects));
binucleateIndexes = find(arrayfun(@(x) x.NumNuclei > 1, Options.Objects));
for n = 1:length(binucleateIndexes) % Label them as excluded.
    Options.Objects(binucleateIndexes(n)).Retained = 0;
end

% Find crops with no objects
noObjectsList = Options.Crop(arrayfun(@(x) (isempty(x.NumCrops) || x.NumCrops == 0), Options.Crop));
for ni = 1:length(noObjectsList) % Label them as excluded.
    objectIndex = getObjectIndex(noObjectsList(ni).File, noObjectsList(ni).Site, 0, Options.Objects);
    for o = objectIndex
        Options.Objects(o).Retained = 0;
    end
end

% Discard rejected if any
if isfield(Options.Objects, 'Rejected')
    for p = 1:length(Options.Objects)
        if isempty(Options.Objects(p).Rejected)
            Options.Objects(p).Rejected = 0;
        elseif Options.Objects(p).Rejected == 1
            Options.Objects(p).Retained = 0;
        elseif Options.Objects(p).Rejected == -1
            Options.Objects(p).Retained = 1;
        end
    end
end

% Remove manually cropped objects that were not segmented because of min
% cell size or min intensity
selectedRows = [];
for p = 1:length(Options.Objects)
    if isempty(Options.Objects(p).NumObjects3D)
        selectedRows = [selectedRows, getObjectIndex(Options.Objects(p).File, Options.Objects(p).Site, Options.Objects(p).Crop, Options.Objects)];
    end
end
Options.Objects(selectedRows) = [];
save([ExpSettings.CloudPath, ExpSettings.ExpID, '-Options.mat'], 'Options', '-v7.3');

for p = 1:length(Objects)
    if isempty(Objects(p).File)
        Objects(p) = [];
    end
end

fileKey = getFileKey(1,0,0,metadata);
save([ExpSettings.ObjectsPath, fileKey, '-Objects.mat'], 'Objects', '-v7.3');
%% Global excluded and retained 

excludedObjectsList = Options.Objects(arrayfun(@(x) x.Retained == 0, Options.Objects));
excludedObjectsData = Objects(arrayfun(@(x) x.Retained == 0, Options.Objects));

retainedObjectsList = Options.Objects(arrayfun(@(x) x.Retained == 1, Options.Objects));
retainedObjectsData = Objects(arrayfun(@(x) x.Retained == 1, Options.Objects));

%%

selectedObjectsList = retainedObjectsList(arrayfun(@(x) x.NumNuclei == 0, retainedObjectsList));
selectedObjectsData = retainedObjectsData(arrayfun(@(x) x.NumNuclei == 0, retainedObjectsList));

% inspectObjects(selectedObjectsList, Options, ExpSettings, metadata);