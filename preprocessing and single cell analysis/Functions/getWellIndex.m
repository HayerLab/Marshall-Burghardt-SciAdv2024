function well = getWellIndex(file, site, metadata)
% Rodrigo Migueles April 2023
    if metadata(file).n_wells > 1
        well = find(site <= cumsum(metadata(file).wells(:,3)), 1);
    else
        well =1;
    end
end