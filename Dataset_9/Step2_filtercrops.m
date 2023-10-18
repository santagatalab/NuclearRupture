function Step2_filtercrops(filename, options) 
    % Constructing the path to the folder containing cropped images
    crop_folder = [filename.datafolder filename.folders.fullstacks];
    
    % Listing the contents of the cropped image folder
    list_crops = dir(crop_folder);
    
    % Loop through each item in the list (starting from index 3 to skip '.' and '..')
    for i = 3:length(list_crops)
        cropname = [crop_folder list_crops(i).name]; % Current cropped image's full path
        
        % Checking if the file extension is either 'tif' or 'tiff'
        if strcmp(cropname(end-2:end), 'tif') || strcmp(cropname(end-3:end), 'tiff')
            % Reading the first index (image) of the current cropped image
            DAPIcheck = imread(cropname, 'Index', 1);
            
            % Checking if the calculated percentile of DAPI intensity is below a threshold
            if prctile(DAPIcheck(:), options.crops.prctile) < options.crops.DAPIbkgd_crops
                % If the condition is met, delete the current cropped image
                delete(cropname)
            end
        end
    end
end