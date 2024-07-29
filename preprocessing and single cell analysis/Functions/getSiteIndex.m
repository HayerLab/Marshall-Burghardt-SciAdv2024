function siteIndex = getSiteIndex(file, site, metadata)
    sitesList = [];
    siteIndex = 0; 
    for f = 1:length(metadata)
        for s = 1:metadata(f).n_sites
            siteIndex = siteIndex + 1;
            sitesList(siteIndex, 1) = siteIndex; %#ok<AGROW>
            sitesList(siteIndex, 2) = f; %#ok<AGROW>
            sitesList(siteIndex, 3) = s; %#ok<AGROW>                        
        end
    end

    siteIndex = (find(sitesList(:,2) == file));
    sitesList = sitesList(siteIndex, :);
    siteIndex = (find(sitesList(:,3) == site));
    sitesList = sitesList(siteIndex, :);
    siteIndex = (find(sitesList(:,3) == site));
    siteIndex = sitesList(siteIndex,1);
end