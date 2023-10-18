function Dataset_2_Step4_CycIF_measurements_TMA_micro(filename, options)  

addpath(filename.codesfolder)
mkdir(filename.resultfolder); 
tic
fullstackfolder = [filename.datafolder 'dearray' filesep]; 
ilastikfolder = [filename.analfolder filename.ilastiksegfol];
ilastikprobfolder = [ilastikfolder 'pmap\'];
ilastiksegfolder = [ilastikfolder 'seg\'];
ilastiklabelfolder = [ilastikfolder 'label\'];
ilastikcheckfolder = [ilastikfolder 'check\'];
ilastiknuccytsegfolder = [ilastikfolder 'nuccytseg\'];
mkdir(ilastiknuccytsegfolder);
if exist([filename.resultfolder filename.resultfile],'file') == 2
    load([filename.resultfolder filename.resultfile])
else
    Results{1} = [];
    Results{1}.CoreFlag = [];
    Results{1}.Area = [];
    Results{1}.Solidity = [];
    Results{1}.CentroidX = [];
    Results{1}.CentroidY = [];
    Results{1}.MedianNucSign = [];
    Results{1}.MedianCytSign = [];
    Results{1}.MeanNucSign = [];
    Results{1}.MeanCytSign = [];
    Results{1}.Name = [];
end

for k1 = 1:length(filename.realcoreinfo)
    rootname = num2str(k1);
    disp([num2str(k1/length(filename.realcoreinfo) *100) ' percent done']);
    if ~isempty(rootname)
        Results{k1} = [];
        Results{k1}.CoreFlag = [];
        Results{k1}.Area = [];
        Results{k1}.Solidity = [];
        Results{k1}.CentroidX = [];
        Results{k1}.CentroidY = [];
        Results{k1}.MedianNucSign = [];
        Results{k1}.MedianCytSign = [];
        Results{k1}.MeanNucSign = [];
        Results{k1}.MeanCytSign = [];
        Results{k1}.Name = rootname;
        %%%% load the FullStack and/to check it exists ----------------
        FileTif = [fullstackfolder rootname  '.ome.tif'];

        if exist(FileTif,'file')~=2
            disp(FileTif)
            disp('Error: FullStack not found')
            continue
        end
        disp(FileTif)

        %%%% load nuclear mask, create cytoplasmic mask and save them -
        %rootname = Results{k1}.Name;
        segfile = [ilastiksegfolder rootname  filename.segsuffix];
        cyt_segfile = [ilastiknuccytsegfolder rootname  '_NucCytSeg.tif'];

        if exist(segfile,'file') ~= 2
            disp(segfile)
            disp('Segmentation file not found')
            continue
        end

        try
            if Results{k1}.CoreFlag == 1
                    disp('Core already analysed')
                    continue
            end
        catch
            disp('Starting Core analysis')
        end
        
        NucMask = uint16(imread(segfile));
        lb_NucMask = uint16(bwlabel(NucMask));

        % create cytoplasmic mask
        lb_CytMask = imdilate(lb_NucMask,offsetstrel('ball',ceil(options.cytosize),0));
        lb_FullCellMask = lb_CytMask;
        lb_CytMask(lb_NucMask>0)=0;
        CytMask = uint16(lb_CytMask > 0);


        DAPI_img = uint16(imread(FileTif,'Index',1));
        NucCytSeg = cat(3,NucMask,DAPI_img,CytMask);
        imwrite(NucCytSeg,cyt_segfile);

        %%%% measure cell properties and assign them to results -------
        stats_nuc = regionprops(lb_NucMask,'Area','Centroid','Solidity');
        stats_cyt = regionprops(lb_CytMask,'Area');

        Area = {stats_nuc.Area};
        Solidity = {stats_nuc.Solidity};
        Centroid = {stats_nuc.Centroid};
        CentroidVect = cell2mat(Centroid);
        CentroidMat = reshape(CentroidVect,2,length(CentroidVect)/2);
        CytArea = {stats_cyt.Area};

        totcells_field = length(Area);
        currentcells = 0; %length(Results{k1}.Area);
        startcell = currentcells + 1;
        endcell = currentcells + totcells_field;

        if totcells_field < 1
            disp(['No cells in Field ' rootname])
            continue
        end

        Results{k1}.Area(startcell:endcell,1) = cell2mat(Area)';
        Results{k1}.Solidity(startcell:endcell,1) = cell2mat(Solidity);
        Results{k1}.CentroidX(startcell:endcell,1) = CentroidMat(1,:);
        Results{k1}.CentroidY(startcell:endcell,1) = CentroidMat(2,:);
        Results{k1}.CytArea(startcell:endcell,1) = [cell2mat(CytArea)'; zeros(totcells_field-length(cell2mat(CytArea)),1)];
           % segment foci or subcellular structures
                   % segment foci or subcellular structures
        if options.substructures.flag == 1 
           i3 = options.substructures.channel;
             FluorImage = imread(FileTif,'Index',i3);
             edgefluor = FluorImage == 0;
             edgeeliminatemask = stdfilt(edgefluor) == 0;
            FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
           Focisv_BK = FluorImage - imopen(FluorImage, strel('disk',options.substructures.bkgsize)); 
           SpatVar_Raw = stdfilt(Focisv_BK)./mean(FluorImage( FluorImage(:)>10  ));

           % FUNDAMENTAL PARAMETER HERE!!!!
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            FociSeg = SpatVar_Raw>options.substructures.SV_thresh;       % we increase this parameter to be more conservative (from 0.8 to 1)
%                 figure(1)
%                 imshow(SpatVar_Raw,[])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % THIS VERSION OF THE ANALYSIS HAS NO DILATION!!!!
            FociSeg = FociSeg > 0 & edgeeliminatemask;
            FociSeg = imfill(FociSeg,'holes');
            FociSeg = bwareaopen(FociSeg,16);
            % THIS VERSION OF THE ANALYSIS HAS NO DILATION!!!!
            Dilate_Nuc_Image = imdilate(lb_NucMask,offsetstrel('ball',options.substructures.nucdilation,0));
            lb_FociImage = uint16(FociSeg).*Dilate_Nuc_Image;
            lb_FociImage = imfill(lb_FociImage,'holes');
            mask = bwareaopen(lb_FociImage,5);
            lb_FociImage = lb_FociImage.*uint16(mask);
            %eccentricity 
            FociSegImage = lb_FociImage > 0;
            Label = labelmatrix(bwconncomp(FociSegImage));
            Nucleiecc = regionprops(FociSegImage,"Circularity");
            thinnuclei = find([Nucleiecc.Circularity] < 0.5);
            for i = thinnuclei
                Label(Label==i) = 0;
            end
            lb_FociImage = lb_FociImage.*uint16(Label>0);
            tempmask = bwareaopen(( (lb_FociImage>0) | (NucMask>0) ),100)>0;
            lb_FociImage = lb_FociImage.*uint16(tempmask);
            FociSegImage = lb_FociImage > 0;

            FociCheckImage = uint16(NucMask);
            FociCheckImage(FociSegImage>0)=2;

            DirectCheckIm = cat(3,uint16(FociCheckImage),uint16(Focisv_BK));
            DirectCheckIm = cat(3,DirectCheckIm,uint16(DAPI_img));


            substructuresegfile = [filename.analfolder filename.ilastiksegfol 'checknucsubseg\' num2str(filename.realcoreinfo(k1).index) '_Seg_SubNucStruct_SV' num2str(options.substructures.SV_thresh) '.tif'];
            imwrite(uint16(DirectCheckIm),substructuresegfile)

            Foci_stats = regionprops(lb_FociImage,FluorImage,'PixelValues');
            Foci_tot_stats = regionprops(Dilate_Nuc_Image,FluorImage,'PixelValues');

            Results{k1}.AreaSubstruct(1:max(lb_FociImage(:)),1) = cellfun(@length,{Foci_stats.PixelValues});
            Results{k1}.AreaSubstruct(1:max(Dilate_Nuc_Image(:)),2) = cellfun(@length,{Foci_tot_stats.PixelValues});

            Results{k1}.MeanInsideSubstruct(1:max(lb_FociImage(:)),1) = cellfun(@sum,{Foci_stats.PixelValues});
            Results{k1}.MeanInsideSubstruct(1:max(Dilate_Nuc_Image(:)),2) = cellfun(@sum,{Foci_tot_stats.PixelValues});

            Results{k1}.AreaSubstruct(isnan(Results{k1}.AreaSubstruct))=0;
            Results{k1}.MeanInsideSubstruct(isnan(Results{k1}.MeanInsideSubstruct))=0;


        end
        %%%% cycle around the images and extract measurements ---------
        for i3 = 1:max(filename.cycles)*4
            % load image
            FluorImage = imread(FileTif,'Index',i3);
            FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
            % extract pixel values from fluo image using the segmented images    
            Nuclei_stats = regionprops(lb_NucMask,FluorImage_BK,'PixelValues');
            Cytopl_stats = regionprops(lb_CytMask,FluorImage_BK,'PixelValues');
            FullCell_stats = regionprops(lb_FullCellMask,FluorImage_BK,'PixelValues');

            Results{k1}.MedianNucSign(startcell:endcell,i3) = cellfun(@median,{Nuclei_stats.PixelValues});
            Results{k1}.MedianCytSign(startcell:endcell,i3) = cellfun(@median,{Cytopl_stats.PixelValues});
            Results{k1}.MeanNucSign(startcell:endcell,i3) = cellfun(@mean,{Nuclei_stats.PixelValues});
            Results{k1}.MeanCytSign(startcell:endcell,i3) = cellfun(@mean,{Cytopl_stats.PixelValues});
            Results{k1}.MeanFullCellSign(startcell:endcell,i3) = cellfun(@mean,{FullCell_stats.PixelValues});
            Results{k1}.MedianFullCellSign(startcell:endcell,i3) = cellfun(@median,{FullCell_stats.PixelValues});

             %All channel Foci segmentation
             Foci_stats = regionprops(lb_FociImage,FluorImage_BK,'PixelValues');
            Results{k1}.MeanFociSign(1:max(lb_FociImage(:)),i3) = cellfun(@mean,{Foci_stats.PixelValues});
            Results{k1}.MedianFociSign(1:max(lb_FociImage(:)),i3) = cellfun(@median,{Foci_stats.PixelValues});
            Results{k1}.MeanFociSign(isnan(Results{k1}.MeanFociSign))=0;
            Results{k1}.MedianFociSign(isnan(Results{k1}.MedianFociSign))=0;
        end
        Results{k1}.CoreFlag = 1; 
        %save([filename.resultfolder filename.resultfile],'Results', '-v7.3') %Saving matrix
        save([filename.resultfolder filename.resultfile],'Results') %Saving matrix
    end
end
