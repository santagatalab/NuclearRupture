function Step1_fields_preilastik(filename, options) 
    % Setting the output folder to the data folder from the filename structure
    outputfold = filename.datafolder; 
    
    % Adding the main folder to the MATLAB path
    filename.folders.main = filename.rawfolder;
    addpath(filename.folders.main)
    
    try
        % Creating the main output folder if it doesn't exist
        mkdir(outputfold) 
        
        % Creating primaryNucfolder and its associated ilastiksegfol subfolder
        mkdir(filename.primaryNucfolder)
        mkdir([filename.primaryNucfolder filename.ilastiksegfol])
        
        % Creating microNucfolder and its associated ilastiksegfol subfolder
        mkdir(filename.microNucfolder)
        mkdir([filename.microNucfolder filename.ilastiksegfol])
    end

    %% Splitting into size of field(t) desired 

    t = filename.sizefield;
    tic
    Coordinates = {}; % Cell array to store field coordinates
    fold = filename.tissue; % Current tissue name
    disp(fold) 

    % Creating folders for the current tissue's outputs
    outputfolder = [outputfold];
    addpath(outputfold)
    mkdir([outputfolder filename.folders.coordinates])  
    mkdir([outputfolder filename.ilastiksegfol])
    mkdir([outputfolder filename.folders.fullstacks])
    
    % Constructing the output filename for saving coordinates
    filenamesv = [outputfold filename.folders.coordinates fold];

    % Constructing paths and filenames for input ometiff file
    omefold = [filename.folders.main];
    omefile = [omefold fold '.ome.tif'];
    
    % Reading the first image of the ometiff file (DAPICycle0)
    DAPICycle0 = imread(omefile, 1);
    [y, x] = size(DAPICycle0); % y = # of pixels in 1 column; x = # of pixels in 1 row
    
    fld_rows = ceil(y / t);
    fld_cols = ceil(x / t);
    
    for r = 1:fld_rows
        for c = 1:fld_cols
            % Constructing field name and output filename
            fieldname = [fold '_Field_' num2str(r, filename.dim) '_' num2str(c, filename.dim)];
            outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
            
            % Defining pixel ranges for current field
            pix_rows = 1 + t * (r - 1):min(y, r * t);
            pix_cols = 1 + t * (c - 1):min(x, c * t);
            
            % Extracting field data and checking a condition
            field = DAPICycle0(pix_rows, pix_cols);
            if prctile(field(:), 97.5) < options.DAPIbkgd
                kf = 0; % Keep field flag
            else
                kf = 1;
                imwrite(field, outfilename)
            end
            % Saving coordinates for the current field
            Coordinates.Field{r, c} = [kf, min(pix_rows), max(pix_rows), min(pix_cols), max(pix_cols)];
        end
    end

    toc
    % Saving the calculated coordinates and other parameters
    save(filenamesv, 'Coordinates', 'x', 'y', 't')
    
    %% Loop through channels to create stacks
    for ch = 2:filename.maxround
        Cycle = imread(omefile, ch);
        
        if size(Cycle) == size(DAPICycle0)
            for r = 1:fld_rows
                for c = 1:fld_cols
                    coordvect = Coordinates.Field{r, c};
                    if coordvect(1) == 1 % Keepfield flag is yes
                        r_coord = coordvect(2:3); % Initial, final
                        c_coord = coordvect(4:5); % Initial, final
                        field = Cycle(r_coord(1):r_coord(2), c_coord(1):c_coord(2));
                        fieldname = [fold '_Field_' num2str(r, filename.dim) '_' num2str(c, filename.dim)];
                        
                        outfilename = [outputfolder filename.folders.fullstacks fieldname '.tif'];
                        imwrite(field, outfilename, 'WriteMode', 'append')
                    end  
                end
            end

            fprintf([num2str(ch) ' '])
        else
            break
        end
    end
    
    fprintf('\n')
    disp(['Tissue ' fold ' done'])
    toc
end