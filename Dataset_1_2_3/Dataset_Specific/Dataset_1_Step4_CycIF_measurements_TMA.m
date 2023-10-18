function Dataset_1_Step4_CycIF_measurements_TMA(filename, options)  

addpath(filename.codesfolder)
mkdir(filename.resultfolder); 
tic
fullstackfolder = [filename.analfolder 'FullStacks' filesep]; 

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
    if ~isempty(filename.realcoreinfo{k1})
        Results{k1} = [];
        
        %%%% load the FullStack and/to check it exists ----------------
        FileTif = [fullstackfolder 'Core'      filename.realcoreinfo{k1}.name  filename.suffix];

        if exist(FileTif,'file')~=2
            disp(FileTif)
            disp('Error: FullStack not found')
            continue
        end
        disp(FileTif)

        %%%% load nuclear mask, create cytoplasmic mask and save them -
        segfile = [filename.analfolder filename.ilastiksegfol 'Core'  filename.realcoreinfo{k1}.name filename.segsuffix];
        cyt_segfile = [filename.analfolder filename.ilastiksegfol 'Core'  filename.realcoreinfo{k1}.name '_NucCytSeg.tif'];

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
        Results{k1}.CoreCoord(startcell:endcell,1) = zeros(totcells_field,1)+find(cell2mat(filename.prefix1) == filename.realcoreinfo{k1}.name(1))*options.tilesize*options.rows;
        Results{k1}.CoreCoord(startcell:endcell,2) = zeros(totcells_field,1)+str2num(filename.realcoreinfo{k1}.name(5:6))*options.tilesize*options.cols;
           % segment foci or subcellular structures
        if options.substructures.flag == 1 
           i3 = options.substructures.channel;
           FluorImage = imread(FileTif,'Index',i3);
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
            FociSeg = FociSeg > 0;
            % THIS VERSION OF THE ANALYSIS HAS NO DILATION!!!!
            Dilate_Nuc_Image = imdilate(lb_NucMask,offsetstrel('ball',options.substructures.nucdilation,0));
            lb_FociImage = uint16(FociSeg).*Dilate_Nuc_Image;
            lb_FociImage = imfill(lb_FociImage,'holes');
            FociSegImage = lb_FociImage > 0;
            lb_nucedgefocimask = lb_FociImage;
            lb_nucedgefocimask(lb_nucedgemask>0) = 0;

            FociCheckImage = uint16(NucMask);
            FociCheckImage(FociSegImage>0)=2;

            DirectCheckIm = cat(3,uint16(FociCheckImage),uint16(Focisv_BK));
            DirectCheckIm = cat(3,DirectCheckIm,uint16(DAPI_img));

            substructuresegfile = [filename.analfolder filename.ilastiksegfol 'Core'  filename.realcoreinfo{k1}.name '_Seg_SubNucStruct_SV' num2str(options.substructures.SV_thresh) '.tif'];
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
%         for i3 = 1:max(filename.cycles)*4
%             % load image
%             FluorImage = imread(FileTif,'Index',i3);
%             FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
%             % extract pixel values from fluo image using the segmented images    
%             Nuclei_stats = regionprops(lb_NucMask,FluorImage_BK,'PixelValues');
%             Cytopl_stats = regionprops(lb_CytMask,FluorImage_BK,'PixelValues');
% 
%             Results{k1}.MedianNucSign(startcell:endcell,i3) = cellfun(@median,{Nuclei_stats.PixelValues});
%             Results{k1}.MedianCytSign(startcell:endcell,i3) = cellfun(@median,{Cytopl_stats.PixelValues});
%             Results{k1}.MeanNucSign(startcell:endcell,i3) = cellfun(@mean,{Nuclei_stats.PixelValues});
%             Results{k1}.MeanCytSign(startcell:endcell,i3) = cellfun(@mean,{Cytopl_stats.PixelValues});
% 
%                %All channel Foci segmentation
%              Foci_stats = regionprops(lb_FociImage,FluorImage_BK,'PixelValues');
%             Results{k1}.MeanFociSign(1:max(lb_FociImage(:)),i3) = cellfun(@mean,{Foci_stats.PixelValues});
%             Results{k1}.MedianFociSign(1:max(lb_FociImage(:)),i3) = cellfun(@median,{Foci_stats.PixelValues});
%             Results{k1}.MeanFociSign(isnan(Results{k1}.MeanFociSign))=0;
%             Results{k1}.MedianFociSign(isnan(Results{k1}.MedianFociSign))=0;
%         end

        for i3 = 1:max(filename.cycles)*4
            % load image
            FluorImage = imread(FileTif,'Index',i3);
%             if i3 == 7
%                 % Filtering out autoflorescent signal
%                 segmask = FluorImage < 600;
%                 FluorImage(segmask) = 0; 
%             end
%             if i3 == 
%                 % Filtering out autoflorescent signal
%                 segmask = FluorImage < 600;
%                 FluorImage(segmask) = 0; 
%             end
%             
            FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
            % extract pixel values from fluo image using the segmented images    
            Nuclei_stats = regionprops(lb_NucMask,FluorImage_BK,'PixelValues');
            Cytopl_stats = regionprops(lb_CytMask,FluorImage_BK,'PixelValues');
            Edge_stats = regionprops(lb_nucedgemask,FluorImage_BK,'PixelValues');

            Results{k1}.MedianNucSign(startcell:endcell,i3) = cellfun(@median,{Nuclei_stats.PixelValues});
            Results{k1}.MedianCytSign(startcell:endcell,i3) = cellfun(@median,{Cytopl_stats.PixelValues});
            Results{k1}.MedianEdgeSign(startcell:endcell,i3) = cellfun(@median,{Edge_stats.PixelValues});
            Results{k1}.MeanNucSign(startcell:endcell,i3) = cellfun(@mean,{Nuclei_stats.PixelValues});
            Results{k1}.MeanCytSign(startcell:endcell,i3) = cellfun(@mean,{Cytopl_stats.PixelValues});
            Results{k1}.MeanEdgeSign(startcell:endcell,i3) = cellfun(@mean,{Edge_stats.PixelValues});

               %All channel Foci segmentation
             Foci_stats = regionprops(lb_FociImage,FluorImage_BK,'PixelValues');
            Results{k1}.MeanFociSign(1:max(lb_FociImage(:)),i3) = cellfun(@mean,{Foci_stats.PixelValues});
            Results{k1}.MedianFociSign(1:max(lb_FociImage(:)),i3) = cellfun(@median,{Foci_stats.PixelValues});
            Results{k1}.MeanFociSign(isnan(Results{k1}.MeanFociSign))=0;
            Results{k1}.MedianFociSign(isnan(Results{k1}.MedianFociSign))=0;
        end
        
        
        
        Results{k1}.CoreFlag = 1; 
        save([filename.resultfolder filename.resultfile],'Results', '-v7.3') %Saving matrix
    end
end
