function Stack2TIFF(data, TIFF_path)
% Stack2TIFF
% Rodrigo Migueles Oct 2022

    for frame = 1:size(data,3)
        if frame == 1; WriteMode = "overwrite"; else; WriteMode = "append"; end
        imwrite(uint16(data(:,:,frame)), TIFF_path,...
            "WriteMode", WriteMode, "Compression", "none");
    end

end