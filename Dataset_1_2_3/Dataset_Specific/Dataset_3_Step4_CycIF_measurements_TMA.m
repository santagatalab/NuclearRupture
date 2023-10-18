function Dataset_3_Step4_CycIF_measurements_TMA(filename, options)  

addpath(filename.codesfolder)
mkdir(filename.resultfolder); 
mkdir([filename.analfolder filename.ilastiksegfol 'checknucsubseg\']);
tic
fullstackfolder = [filename.datafolder 'dearray' filesep];
ilastikfolder = [filename.analfolder filename.ilastiksegfol];
ilastiksegfolder = [ilastikfolder 'seg\'];
ilastikNucCytSegfolder = [ilastikfolder 'nuccytseg\'];
mkdir(ilastikNucCytSegfolder);

if exist([filename.resultfolder filename.resultfile],'file') == 2
    load([filename.resultfolder filename.resultfile])
else
    Results{1} = [];
    Results{1}.CoreCoord = [];
    Results{1}.CoreFlag = [];
    Results{1}.Area = [];
    Results{1}.Solidity = [];
    Results{1}.CentroidX = [];
    Results{1}.CentroidY = [];
    Results{1}.MedianNucSign = [];
    Results{1}.MedianCytSign = [];
    Results{1}.MeanNucSign = [];
    Results{1}.MeanCytSign = [];
    if options.substructures.flag == 1
        Results{1}.AreaSubstruct = [];
        Results{1}.MeanInsideSubstruct = [];
%         Results.MeanOutsideSubstruct = [];
%         Results.MedianInsideSubstruct = [];
%         Results.MedianOutsideSubstruct = [];
    end
end
        

for k1 = 1:length(filename.realcoreinfo)
    totallength = length(filename.realcoreinfo);
    disp([num2str(k1/totallength *100) ' percent done']);
    FileTif = [fullstackfolder num2str(filename.realcoreinfo(k1).index)  '.ome.tif'];
    if ~isempty(FileTif)
        if exist(FileTif,'file')~=2
            disp(FileTif)
            disp('Error: FullStack not found')
            continue
        end
        disp(FileTif)

        %%%% load nuclear mask, create cytoplasmic mask and save them -
        segfile = [ilastiksegfolder num2str(filename.realcoreinfo(k1).index) filename.segsuffix];
        cyt_segfile = [ilastikNucCytSegfolder num2str(filename.realcoreinfo(k1).index) '_NucCytSeg.tif'];

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

        lb_innermask = imerode(lb_NucMask,offsetstrel('ball',3,0));
        lb_nucedgemask = lb_NucMask; 
        lb_nucedgemask( lb_innermask > 0) = 0;


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
            disp(['No cells in Field ' FieldFile])
            continue
        end

        Results{k1}.Area(startcell:endcell,1) = cell2mat(Area)';
        Results{k1}.Solidity(startcell:endcell,1) = cell2mat(Solidity);
        Results{k1}.CentroidX(startcell:endcell,1) = CentroidMat(1,:);
        Results{k1}.CentroidY(startcell:endcell,1) = CentroidMat(2,:);
        Results{k1}.CytArea(startcell:endcell,1) = [cell2mat(CytArea)'; zeros(totcells_field-length(cell2mat(CytArea)),1)];
           % segment foci or subcellular structures
        if options.substructures.flag == 1 
           i3 = options.substructures.channel;
             FluorImage = imread(FileTif,'Index',i3);
             edgefluor = FluorImage == 0;
             edgeeliminatemask = stdfilt(edgefluor) == 0;
            FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*1)));
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
            % THIS VERSION OF THE ANALYSIS HAS NO DILATION!!!!
            Dilate_Nuc_Image = imdilate(lb_NucMask,offsetstrel('ball',options.substructures.nucdilation,0));
            %Inner_Nuc_Image = imerode(lb_NucMask,offsetstrel('ball',options.substructures.blockinner,0));
            lb_FociImage = uint16(FociSeg).*Dilate_Nuc_Image;
            %lb_FociImage(Inner_Nuc_Image>0) = -1;
            lb_FociImage = imfill(lb_FociImage,'holes');
            %micronucimage = imread(['Z:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2022_11_GBM399_TMA_Lineage\LSP14356_Rarecyte\ANALYSIS20230130_micro\seg\' rootname '_Seg.tif']);
            %micromask = micronucimage > 0;
            %micromask = imdilate(micromask,strel('square',3));
            %lb_FociImage(micromask) = -1;
            mask = bwareaopen(lb_FociImage,10);
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

            %remove apoptotic cell by area
            Fociarea = regionprops(lb_FociImage,"Area");
            Fociarea = [Fociarea.Area];
            Nucleararea = regionprops(Dilate_Nuc_Image,"Area");
            Nucleararea = [Nucleararea.Area];
            Fociarea(length(Fociarea)+1:length(Nucleararea)) = 0;
            Arearatio = Fociarea > (0.5 * Nucleararea);
            Arearatioindex = find(Arearatio);
            %FociSegImage = lb_FociImage > 0;
            %Label = labelmatrix(bwconncomp(FociSegImage));
            for i = Arearatioindex
                %Label(Label==i) = 0;
                lb_FociImage(lb_FociImage==i) = 0;
            end
            %lb_FociImage = lb_FociImage.*uint16(Label>0);
            FociSegImage = lb_FociImage > 0;

            FociCheckImage = uint16(NucMask);
            FociCheckImage(FociSegImage>0)=2;

            DirectCheckIm = cat(3,uint16(FociCheckImage),uint16(Focisv_BK));
            DirectCheckIm = cat(3,DirectCheckIm,uint16(DAPI_img));

            substructuresegfile = [ilastikfolder 'checknucsubseg\' num2str(filename.realcoreinfo(k1).index) '_Seg_SubNucStruct_SV' num2str(options.substructures.SV_thresh) '.tif'];
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
            %if ~ismember(ceil(i3/4), filename.cycles)
            %    continue
            %end
            %if i3 == 2
            %    continue
            %end
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
            Results{k1}.MedianFociSign( size(Results{k1}.Area,1),1) = uint16(0);
            Results{k1}.MeanFociSign( size(Results{k1}.Area,1),1) = 0;
        end
        Results{k1}.CoreFlag = 1; 
        %save([filename.resultfolder filename.resultfile],'Results', '-v7.3') %Saving matrix
        save([filename.resultfolder filename.resultfile],'Results') %Saving matrix
    end
end
for i = 1:length(Results)
    if ~isempty(Results{i})
        Results{i}.Corenum(1:size(Results{i}.Area,1),1) = i;
        if size(Results{i}.MeanFociSign,1) < size(Results{i}.Area,1) %Correction for NaN being deleted 
            Results{i}.MedianFociSign( size(Results{i}.Area,1),1) = uint16(0);
            Results{i}.MeanFociSign( size(Results{i}.Area,1),1) = 0;
        end
    end
end
 save([filename.resultfolder filename.resultfile],'Results', '-v7.3') %Saving matrix