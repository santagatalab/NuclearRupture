function Step3_segmentfromilastik_micro(filename, options) 
    % Open the tiff file of the full stack images and segments them based on Ilastik Probabilities file
    % Inputs: FullStacks and Ilastik Probabilities files for micronucleus
    tic
    fullstackfolder = [filename.datafolder 'FullStacks' filesep]; 
    ilastikfolder = [filename.microNucfolder filename.ilastiksegfol];
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
            IlastikTif = [ilastikprobfolder filename.realcoreinfo(i1).name(1:end-4) filename.ilastiksuffix ];
            segfile = [ilastiksegfolder num2str(filename.realcoreinfo(i1).index) filename.segsuffix];
            labelfile = [ilastiklabelfolder num2str(filename.realcoreinfo(i1).index) '_label.tif'];
            
            % Read necessary images and data
            DAPI = imread(FileTif,'Index',1);
            DAPI_last = imread(FileTif,'Index',max(filename.cycles)*4-3);
            IlastikRGB = imread(IlastikTif);
            
            NucProb = IlastikRGB(:,:,options.nuc);
            CytProb = IlastikRGB(:,:,options.cyt);
            BkgProb = IlastikRGB(:,:,options.backgr);
    
            % Pre-filtering based on probability thresholds
            nuc_mask = zeros(size(NucProb));
            temp_mask = zeros(size(NucProb)); 
            else_mask = zeros(size(NucProb));
    
            index = BkgProb > options.max_prob*options.bkgdprob_min | CytProb > options.max_prob*options.cytprob_min; 
            else_mask(index) = 1;
    
            index = NucProb > options.max_prob*options.nucprob_min_micro; 
            temp_mask(index) = 1;
    
            nuc_mask = temp_mask & ~else_mask;
            nuc_mask = bwareaopen(nuc_mask, 5, 4);
    
            sImage = NucProb;
            ThresImage = nuc_mask;
            
            % Import primary nucleus segmentation
            primaryseg = imread([filename.primaryNucfolder 'IlastikSegmentation\seg\' num2str(i1) '_Seg.tif']);
            primarymask = primaryseg > 0;
            joinedmask = ThresImage | primarymask;
            largemask = bwareaopen(joinedmask, 180, 4);
            ThresImage(largemask) = 0;
            largemask = bwareaopen(ThresImage, 80, 4);
            ThresImage(largemask) = 0;
    
            ThresImage = imclose(ThresImage, strel('disk', 1));
            ThresImage = imfill(ThresImage, 'holes');
    
            % Remove thin nuclei based on circularity
            Label = labelmatrix(bwconncomp(ThresImage));
            Nucleiecc = regionprops(ThresImage, "Circularity");
            thinnuclei = find([Nucleiecc.Circularity] < 1);
            for i = thinnuclei
                Label(Label==i) = 0;
            end
            ThresImage = ThresImage & Label;
    
            % Remove false positives with faint background 
            Focisv_BK = DAPI - imopen(DAPI, strel('disk', options.substructures.bkgsize)); 
            SpatVar_Raw = stdfilt(Focisv_BK) ./ mean(DAPI( DAPI(:) > 10 ));
            Label = labelmatrix(bwconncomp(ThresImage));
            for i = 1:max(max(Label))
                mask = Label == i;
                dilutemask = imdilate(mask, strel('disk', 1));
                outmask = dilutemask & ~mask;
                meaninside = mean(DAPI(mask));
                meanoutside = mean(DAPI(outmask));
                if meaninside < 1.2 * meanoutside || meaninside < 2000 || prctile(SpatVar_Raw(mask | outmask), 95) < 0.35
                    Label(Label == i) = 0;
                end
            end
            ThresImage = ThresImage & Label;
            
            % Remove nuclei touching edges
            [mx, my] = size(ThresImage);
            Label = labelmatrix(bwconncomp(ThresImage));
            rp = regionprops(ThresImage);
            ct = {rp.Centroid};
            edge = 20;
            for i = 1:length(ct)
                cx = ct{i}(1);
                cy = ct{i}(2);
                if cx < edge || cx > (mx - edge) || cy < edge || cy > (my - edge)
                    Label(Label == i) = 0;
                    continue
                end
            end
            ThresImage = ThresImage & Label;
    
            % Save checked nuclei image for verification
            Nuclei = ThresImage;
            check_img = cat(3, uint16(Nuclei), uint16(nuc_mask & ~Nuclei), uint16(DAPI));
            check_file = [ilastikcheckfolder num2str(filename.realcoreinfo(i1).index) '_check' filename.segsuffix];
            imwrite(check_img, check_file)
    
            % Label the nuclei
            imwrite(uint16(Nuclei), segfile)
            label = labelmatrix(bwconncomp(Nuclei));
            imwrite(label, labelfile)
    
        end
    end
end




