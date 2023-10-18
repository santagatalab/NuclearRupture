function Dataset_1_Step3_segmentfromilastik_TMA(filename, options) 
% open the tiff file of the full stack images and segments them based on
% Ilastik Probabilities file 
% Inputs FullStacks and Ilastik Probabilities files 

tic
ilastiksegfol = [filename.analfolder filename.ilastiksegfol];
mkdir(ilastiksegfol)
fullstackfolder = [filename.analfolder 'FullStacks' filesep]; 

for i1 = 1:length(filename.realcoreinfo)
    if ~isempty(filename.realcoreinfo{i1})
        FileTif = [fullstackfolder 'Core'      filename.realcoreinfo{i1}.name  filename.suffix];
        IlastikTif = [fullstackfolder 'Core'      filename.realcoreinfo{i1}.name filename.ilastiksuffix ];
        segfile = [ilastiksegfol 'Core'      filename.realcoreinfo{i1}.name filename.segsuffix];
        labelfile = [ilastiksegfol 'Core'      filename.realcoreinfo{i1}.name '_label'];
        
        if exist(segfile,'file')==2
            disp(['Image ' core ' alredy segmented'])
            continue
        end
            
        disp(FileTif)

        try 
            DAPI = imread(FileTif,'Index',1);
            DAPI_last = imread(FileTif,'Index',max(filename.cycles)*4-3);
        catch
            disp('Error: DAPI not found')
            continue
        end
        

        try
            IlastikRGB = imread(IlastikTif);
        catch
            disp('Error: Ilastik RGB not found')
            disp(IlastikTif)
            continue
        end
        
        NucProb = IlastikRGB(:,:,options.nuc);
        CytProb = IlastikRGB(:,:,options.cyt);
        BkgProb = IlastikRGB(:,:,options.backgr);


        % pre-filter easy ones
        nuc_mask = zeros(size(NucProb));
        temp_mask = zeros(size(NucProb));
        else_mask = zeros(size(NucProb));

        index = BkgProb > options.max_prob*options.bkgdprob_min | CytProb > options.max_prob*options.cytprob_min; 
        else_mask(index) = 1;

        index = NucProb > options.max_prob*options.nucprob_min; 
        temp_mask(index) = 1;

        nuc_mask = temp_mask & ~else_mask;
        nuc_mask = bwareaopen(nuc_mask,options.cellsize*2,4);

        sImage = NucProb;
        ThresImage=nuc_mask;

        seeds = imclose(sImage,strel('disk',2));
        seeds = imgaussfilt(seeds,5,'FilterSize',round(options.cellsize));
        seeds = imregionalmax(seeds);
        seeds=seeds&ThresImage;
        seeds=imdilate(seeds,strel('disk',2));

        % threshold and watershed image
        L=watershed(imimposemin(-sImage,seeds));
        Nuclei=ThresImage&L;

        Nuclei = bwareaopen(Nuclei,options.cellsize*2,4);
        Nuclei = imfill(Nuclei,'holes');
        
        check_img = cat(3,uint16(Nuclei),uint16(DAPI),uint16(DAPI_last));
        check_file = [ilastiksegfol 'Core'      filename.realcoreinfo{i1}.name  '_check' filename.segsuffix];
        imwrite(check_img,check_file)

        % label the nuclei

        imwrite(uint16(Nuclei),segfile)
        imwrite(L,labelfile)

        if options.pausetime > 0
            disp(['Pausing for ' num2str(options.pausetime) ' seconds'])
            pause(options.pausetime)
        end
        
    end
end

