function Dataset_3_Step3_segmentfromilastik_TMA(filename, options) 
% open the tiff file of the full stack images and segments them based on
% Ilastik Probabilities file 
% Inputs FullStacks and Ilastik Probabilities files 

tic
fullstackfolder = filename.fullstackfolder;
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
    FileTif = [fullstackfolder filename.realcoreinfo(i1).name];
    if exist(FileTif,'file') 
        IlastikTif = [ilastikprobfolder filename.realcoreinfo(i1).name(1:end-8) filename.ilastiksuffix ];
        segfile = [ilastiksegfolder num2str(filename.realcoreinfo(i1).index) filename.segsuffix];
        labelfile = [ilastiklabelfolder num2str(filename.realcoreinfo(i1).index) '_label.tif'];

        try 
            DAPI = imread(FileTif,'Index',1);
            DAPI_last = imread(FileTif,'Index',max(filename.cycles)*4-3);
        catch
            disp('Error: DAPI not found')
%             continue
        end
        

        try
            IlastikRGB = imread(IlastikTif,1);
        catch
            disp('Error: Ilastik RGB not found')
            disp(IlastikTif)
            continue
        end
        
        NucProb = IlastikRGB(:,:,options.nuc);
        CytProb = IlastikRGB(:,:,options.cyt);
        BkgProb = IlastikRGB(:,:,options.backgr);

        %NucProb = imopen(NucProb,strel('octagon',3));
        %BkgProb = imopen(BkgProb,strel('octagon',3));
        
        % pre-filter easy ones
        nuc_mask = zeros(size(NucProb));
        temp_mask = zeros(size(NucProb));
        else_mask = zeros(size(NucProb));

        index = BkgProb > options.max_prob*options.bkgdprob_min | CytProb > options.max_prob*options.cytprob_min; 
        else_mask(index) = 1;

        index = NucProb > options.max_prob*options.nucprob_min; 
        temp_mask(index) = 1;

        %nuc_mask = temp_mask & ~else_mask;
        %nuc_mask = bwareaopen(nuc_mask,options.cellsize*2,4);

        nuc_mask = temp_mask & ~else_mask;
        nuc_mask = bwareaopen(nuc_mask,100);
        %nuc_mask = imfill(nuc_mask,'holes');
        nuc_mask = imdilate(nuc_mask,strel('disk', 3));
        
        

        sImage = NucProb;
        ThresImage=nuc_mask;

        seeds = imclose(sImage,strel('disk',2));
        seeds = imgaussfilt(seeds,5,'FilterSize',round(options.cellsize/2*1)*2+1);
        seeds = imregionalmax(seeds);
        seeds=seeds&ThresImage;
        seeds=imdilate(seeds,strel('disk',3));

        % threshold and watershed image
        blurry = imclose(sImage,strel('disk',3));
        blurry = imgaussfilt(blurry,5,'FilterSize',15);
        L=watershed(imimposemin(imcomplement(blurry),~ThresImage | seeds));
        Nuclei=ThresImage&L;

        %L=watershed(imimposemin(-sImage,seeds));
        %Nuclei=uint32(ThresImage&L);
        
        %Nuclei = imfill(Nuclei,'holes');
        %Nuclei = bwareaopen(Nuclei,175,4);
        %Nuclei = bwareaopen(Nuclei,options.cellsize*2,4);
        Nuclei = bwareaopen(Nuclei,110,4);
        Nuclei = imfill(Nuclei,'holes');
        forcheck = uint16(Nuclei);
        checkNuclei = bwlabel(uint32(Nuclei));
%         
%         Nuclei = imdilate(Nuclei,strel('disk',2));
%         checkNuclei = imdilate(checkNuclei,strel('disk',2));
        
        
        %check_img = cat(3,uint16(Nuclei),uint16(DAPI),uint16(DAPI_last));
        check_img = cat(3,uint16(checkNuclei),uint16(forcheck),uint16(DAPI));
        %temp_img = cat(3,uint16(else_mask),uint16(DAPI),uint16(forcheck));
        check_file = [ilastikcheckfolder num2str(filename.realcoreinfo(i1).index) '_check' filename.segsuffix];
        %temp_check_file = [ilastikcheckfolder num2str(filename.realcoreinfo(i1).index) '_tempcheck' filename.segsuffix];
        imwrite(check_img,check_file)
        %imwrite(temp_img,temp_check_file)

        % label the nuclei

        %saveastiff(uint32(Nuclei),segfile)
        imwrite(uint16(Nuclei),segfile)
        imwrite(L,labelfile)
        
        
    end
end
% 
% %testing
% im = cat(3,Nuclei,seedsave,Nuclei);
% a= zeros(size(im));
% a(im) = 60000;
% imshow(a)