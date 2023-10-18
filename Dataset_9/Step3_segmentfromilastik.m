function Step3_segmentfromilastik(filename, options) 
    % open the tiff file of the full stack images and segments them based on
    % Ilastik Probabilities file 
    % Inputs FullStacks and Ilastik Probabilities files 
    fullstackfolder = filename.fullstackfolder;
    ilastikfolder = [filename.primaryNucfolder filename.ilastiksegfol];
    mkdir(ilastikfolder)
    ilastikprobfolder = [ilastikfolder 'pmap\'];
    mkdir(ilastikprobfolder)
    ilastiksegfolder = [ilastikfolder 'seg\'];
    mkdir(ilastiksegfolder)
    ilastiklabelfolder = [ilastikfolder 'label\'];
    mkdir(ilastiklabelfolder)
    ilastikcheckfolder = [ilastikfolder 'check\'];
    mkdir(ilastikcheckfolder)
    
    % Loop through each core's information
    for i1 = 1:length(filename.realcoreinfo)
        FileTif = [fullstackfolder filename.realcoreinfo(i1).name];
        
        % Check if the core file exists
        if exist(FileTif,'file') 
            
            % Constructing paths for output files
            IlastikTif = [ilastikprobfolder filename.realcoreinfo(i1).name(1:end-4) filename.ilastiksuffix];
            segfile = [ilastiksegfolder num2str(filename.realcoreinfo(i1).index) filename.segsuffix];
            labelfile = [ilastiklabelfolder num2str(filename.realcoreinfo(i1).index) '_label.tif'];
            
            % Check if the segmentation file already exists
            if exist(segfile,'file')==2
                disp(['Image ' core ' already segmented'])
                continue
            end
            
            disp(FileTif)
    
            try 
                DAPI = imread(FileTif,'Index',1);
                DAPI_last = imread(FileTif,'Index',max(filename.cycles)*4-3);
            catch
                disp('Error: DAPI not found')
                % continue
            end
            
            try
                IlastikRGB = imread(IlastikTif,1);
            catch
                disp('Error: Ilastik RGB not found')
                disp(IlastikTif)
                continue
            end
            
            %Seperating Probability files
            NucProb = IlastikRGB(:,:,options.nuc);
            CytProb = IlastikRGB(:,:,options.cyt);
            BkgProb = IlastikRGB(:,:,options.backgr);
            
            % Creating masks based on probability thresholds
            nuc_mask = zeros(size(NucProb));
            temp_mask = zeros(size(NucProb));
            else_mask = zeros(size(NucProb));
            index = BkgProb > options.max_prob * options.bkgdprob_min | CytProb > options.max_prob * options.cytprob_min; 
            else_mask(index) = 1;
            index = NucProb > options.max_prob * options.nucprob_min; 
            temp_mask(index) = 1;
            nuc_mask = temp_mask & ~else_mask;
            nuc_mask = bwareaopen(nuc_mask, options.cellsize+10);
            nuc_mask = imdilate(nuc_mask, strel('disk', 3));
 
            % Threshold and watershed the image
            sImage = NucProb;
            ThresImage = nuc_mask;
            seeds = imclose(sImage, strel('disk', 2));
            seeds = imgaussfilt(seeds, 5, 'FilterSize', round(options.cellsize/4 ) * 2 + 1);
            seeds = imregionalmax(seeds);
            seeds = seeds & ThresImage;
            seeds = imdilate(seeds, strel('disk', 3));
            blurry = imclose(sImage, strel('disk', 3));
            blurry = imgaussfilt(blurry, 5, 'FilterSize', 15);
            L = watershed(imimposemin(imcomplement(blurry), ~ThresImage | seeds));
            Nuclei = ThresImage & L;
            Nuclei = bwareaopen(Nuclei, options.cellsize-10);

            % Save checked nuclei image for verification
            forcheck = uint16(Nuclei);
            checkNuclei = bwlabel(uint32(Nuclei));
            check_img = cat(3, uint16(checkNuclei), uint16(forcheck), uint16(DAPI));
            check_file = [ilastikcheckfolder num2str(filename.realcoreinfo(i1).index) '_check' filename.segsuffix];
            imwrite(check_img, check_file)
            
            % Save labeled nuclei and segmentation images
            imwrite(uint16(Nuclei), segfile)
            imwrite(L, labelfile)
        end
    end
end
