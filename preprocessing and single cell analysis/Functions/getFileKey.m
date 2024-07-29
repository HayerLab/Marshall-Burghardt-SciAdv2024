function [file_key] = getFileKey(file, site, channel, metadata, convention, well, cropIndex)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
    if nargin < 7
        cropIndex = 0;
        if nargin < 6
            well = getWellIndex(file, site, metadata);
            if nargin < 5
                convention = 'R';
            end       
        end
    end

    if strcmp(convention, 'R')        
        sep = '-';
    elseif strcmp(convention, 'S')
        sep = '_';
    end
    
    site = double(site);

    if metadata(file).n_sites < 100 && site > 0
        siteName = [sep, num2str(site/100),'0']; siteName = [siteName(1), siteName(4:5)];
    elseif metadata(file).n_sites < 999 && metadata(file).n_sites > 99 && site < 100 && site > 9
        siteName = [sep, '0', num2str(site)];
    elseif metadata(file).n_sites < 999 && metadata(file).n_sites > 99 && site < 10 && site > 0
        siteName = [sep, '00', num2str(site)];
    elseif metadata(file).n_sites < 999 && metadata(file).n_sites > 99 && site > 99            
        siteName = [sep, num2str(site)];
    elseif metadata(file).n_sites < 99 && site > 9
        siteName = [sep, num2str(site)];
    elseif metadata(file).n_sites < 9 && metadata(file).n_sites > 0
        siteName = [sep, '0', num2str(site)];            
    elseif site == 0
        siteName = '';
    end

    if cropIndex > 0 && cropIndex < 10
        cropName = [sep, '0', num2str(cropIndex)];
    elseif cropIndex >= 10 && cropIndex < 100
        cropName = [sep, num2str(cropIndex)];            
    else
        cropName = '';
    end
    
    if strcmp(convention, 'R')
        if metadata(file).n_wells > 1 && well > 0
            if channel ~= 0
                file_key = [char(metadata(file).file_name), siteName, cropName, sep, char(metadata(file).wellNames{well, 1}), sep, char(metadata(file).channels{channel, 1})];
            else
                file_key = [char(metadata(file).file_name), siteName, cropName, sep, char(metadata(file).wellNames{well, 1})];
            end        
        else
            if channel ~= 0
                file_key = [char(metadata(file).file_name), siteName, cropName, sep, char(metadata(file).channels{channel, 1})];
            else
                file_key = [char(metadata(file).file_name), siteName, cropName];
            end
        end

    elseif strcmp(convention, 'S')
        if metadata(file).n_rows && metadata(file).n_cols == 1
            row = 1; col = 1; 
            if channel ~= 0
                file_key = [num2str(row), sep, num2str(col), sep, num2str(site), cropName sep, char(metadata(file).channels{channel, 1})];
            else
                file_key = [num2str(row), sep, num2str(col), sep, num2str(site), cropName];
            end            
        else
            [row, col] = getRowColIndex(file, site, metadata);
            if ~isempty(metadata(file).wellNames{well, 1})
                if channel ~= 0
                    file_key = [num2str(row), sep, num2str(col), sep, num2str(site), cropName, sep, char(metadata(file).wellNames{well, 1}), sep, char(metadata(file).channels{channel, 1})];
                else
                    file_key = [num2str(row), sep, num2str(col), sep, num2str(site), cropName, sep, char(metadata(file).wellNames{well, 1})];
                end
            else
                if channel ~= 0
                    file_key = [num2str(row), sep, num2str(col), sep, num2str(site), cropName, sep, char(metadata(file).channels{channel, 1})];
                else
                    file_key = [num2str(row), sep, num2str(col), sep, num2str(site), cropName];
                end                
            end
        end
        
    end

end

