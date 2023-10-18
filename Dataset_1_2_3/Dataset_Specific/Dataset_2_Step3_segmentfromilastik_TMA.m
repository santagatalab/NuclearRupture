function Dataset_2_Step3_segmentfromilastik_TMA(filename, options) 
% open the tiff file of the full stack images and segments them based on
% Ilastik Probabilities file 
% Inputs FullStacks and Ilastik Probabilities files 
%micronucleus
tic
fullstackfolder = [filename.datafolder 'dearray' filesep]; 
ilastikfolder = [filename.analfolder filename.ilastiksegfol];
mkdir(ilastikfolder)
ilastikprobfolder = [ilastikfolder 'pmap\'];
mkdir(ilastikprobfolder)
ilastiksegfolder = [ilastikfolder 'seg\'];
mkdir(ilastiksegfolder)
ilastiklabelfolder = [ilastikfolder 'label\'];
mkdir(ilastiklabelfolder)
ilastikcheckfolder = [ilastikfolder 'check\'];
mkdir(ilastikcheckfolder)

for i1 = 1:length(filename.realcoreinfo)
    FileTif = [fullstackfolder num2str(filename.realcoreinfo(i1).index)  '.ome.tif'];
    if exist(FileTif,'file')
        IlastikTif = [ilastikprobfolder num2str(filename.realcoreinfo(i1).index) filename.ilastiksuffix ];
        segfile = [ilastiksegfolder num2str(filename.realcoreinfo(i1).index) filename.segsuffix];
        labelfile = [ilastiklabelfolder num2str(filename.realcoreinfo(i1).index) '_label.tif'];
        
        DAPI = imread(FileTif,'Index',1);
        DAPI_last = imread(FileTif,'Index',max(filename.cycles)*4-3);
        IlastikRGB = imread(IlastikTif);
        
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
        nuc_mask = bwareaopen(nuc_mask,9,4);

        sImage = NucProb;
        ThresImage=nuc_mask;
        %Import primary nucleus
        primaryseg = imread(['Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Inteferons\ANALYSIS03302023\IlastikSegmentation\seg\' num2str(i1) '_Seg.tif']);
        primarymask = primaryseg>0;
        %primarymask = imdilate(primarymask, strel('disk',3));
        joinedmask = ThresImage | primarymask;
        largemask = bwareaopen(joinedmask,110,4);
        ThresImage(largemask) = 0;
%         largenucleus = NucProb.*uint8(ThresImagemask);
%         largenucleusmask = largenucleus>options.max_prob*options.nucprob_min*1.2; 
%         largenucleusmask = imdilate(largenucleusmask, strel('disk',1));
%         mask1 = bwareaopen(largenucleusmask,800,4);
%         largenucleusmask(mask1) = 0;
%         largenucleusmask = bwareaopen(largenucleusmask,24,4);
%         ThresImage = ThresImage|largenucleusmask;

        ThresImage = imclose(ThresImage,strel('disk',1));
        ThresImage = imfill(ThresImage,'holes');
        

%         % threshold and watershed image
%         L=watershed(imimposemin(imcomplement(sImage),~nuc_mask | seeds));
%         %L=watershed(ThresImage);

        %eccentricity 
        Label = labelmatrix(bwconncomp(ThresImage));
        Nucleiecc = regionprops(ThresImage,"Circularity");
        thinnuclei = find([Nucleiecc.Circularity] < .6);
        for i = thinnuclei
            Label(Label==i) = 0;
        end
        ThresImage=ThresImage&Label;

        %remove false positive with faint background 
        Label = labelmatrix(bwconncomp(ThresImage));
        for i = 1:max(max(Label))
            mask = Label==i;
            dilutemask = imdilate(mask,strel('disk',1));
            outmask = dilutemask & ~mask;
            meaninside = mean(DAPI(mask));
            meanoutside = mean(DAPI(outmask));
            if meaninside < 1.14*meanoutside
                Label(Label==i) = 0;
            end
        end
        ThresImage=ThresImage&Label;

        [mx,my]= size(ThresImage);
        Label = labelmatrix(bwconncomp(ThresImage));
        rp = regionprops(ThresImage);
        ct = {rp.Centroid};
        for i = 1:length(ct)
            cx = ct{i}(1);
            cy = ct{i}(2);
            if cx < 100 || cx > (mx - 50)
                Label(Label==i) = 0;
                continue
            end
            if cy < 100 || cy > (my - 50)
                Label(Label==i) = 0;
                continue
            end
        end
        Nuclei = ThresImage&Label;



        check_img = cat(3,uint16(Nuclei),uint16(nuc_mask & ~Nuclei),uint16(DAPI));
        check_file = [ilastikcheckfolder num2str(filename.realcoreinfo(i1).index) '_check' filename.segsuffix];
        imwrite(check_img,check_file)

        % label the nuclei
        imwrite(uint16(Nuclei),segfile)
        label = labelmatrix(bwconncomp(Nuclei));
        imwrite(label,labelfile)

        if options.pausetime > 0
            disp(['Pausing for ' num2str(options.pausetime) ' seconds'])
            pause(options.pausetime)
        end
        
    end
end



