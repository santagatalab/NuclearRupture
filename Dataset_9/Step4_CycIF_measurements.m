function Step4_CycIF_measurements(filename, options)  
    % Create the result folder if it doesn't exist
    mkdir(filename.resultfolder);
    
    % Create subfolder for nuclear segmentation check
    mkdir([filename.primaryNucfolder filename.ilastiksegfol 'checknucsubseg\']);
    % Assign folder paths
    fullstackfolder = filename.fullstackfolder;
    ilastikfolder = [filename.primaryNucfolder filename.ilastiksegfol];
    ilastiksegfolder = [ilastikfolder 'seg\'];
    ilastikNucCytSegfolder = [ilastikfolder 'nuccytseg\'];
    mkdir(ilastikNucCytSegfolder);
    % Check if result file exists, if yes, load the results, 
    % Otherwise initialize an empty Results cell
    if exist([filename.resultfolder filename.resultfile],'file') == 2
        load([filename.resultfolder filename.resultfile])
    else
        Results{1} = [];
    end
       
    % Loop through each core
    for k1 = 1:length(filename.realcoreinfo)
        totallength = length(filename.realcoreinfo);
        disp([num2str(k1/totallength *100) ' percent done']);
        FileTif = [fullstackfolder filename.realcoreinfo(k1).name];
        % Check if FullStack image exists
        if ~isempty(FileTif)
            if exist(FileTif,'file')~=2
                disp(FileTif)
                disp('Error: FullStack not found')
                continue
            end
            disp(FileTif)
    
            % load nuclear mask, create cytoplasmic mask and save them -
            segfile = [ilastiksegfolder num2str(filename.realcoreinfo(k1).index) filename.segsuffix];
            cyt_segfile = [ilastikNucCytSegfolder num2str(filename.realcoreinfo(k1).index) '_NucCytSeg.tif'];
    
            % Check if Segmentation file exists
            if exist(segfile,'file') ~= 2
                disp(segfile)
                disp('Segmentation file not found')
                continue
            end
    
            % Check if Core is already analyzed
            try
                if Results{k1}.CoreFlag == 1
                        disp('Core already analyzed')
                        continue
                end
            catch
                disp('Starting Core analysis')
            end
            
            % Read segmentation file
            NucMask = uint16(imread(segfile));
            lb_NucMask = uint16(bwlabel(NucMask));
    
            % create cytoplasmic mask
            lb_CytMask = imdilate(lb_NucMask,offsetstrel('ball',ceil(options.cytosize),1));
            lb_FullCellMask = lb_CytMask;
            lb_FullCellMask(lb_NucMask>0) = lb_NucMask(lb_NucMask>0);
            lb_CytMask(lb_NucMask>0)=0;
            CytMask = uint16(lb_CytMask > 0);
    
            % Load DAPI mask
            DAPI_img = uint16(imread(FileTif,'Index',1));
            NucCytSeg = cat(3,NucMask,DAPI_img,CytMask);
            imwrite(NucCytSeg,cyt_segfile);
    
            %%%% measure cell properties and assign them to results -------
            stats_nuc = regionprops(lb_NucMask,'Area','Centroid','Solidity');
            stats_cyt = regionprops(lb_CytMask,'Area');
    
            % Read Nucleus information
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
    
            % Check if there is any nucleus
            if totcells_field < 1
                disp(['No cells in Field '])
                continue
            end
    
            % Add information to Result Data
            Results{k1}.Area(startcell:endcell,1) = cell2mat(Area)';
            Results{k1}.Solidity(startcell:endcell,1) = cell2mat(Solidity);
            Results{k1}.CentroidX(startcell:endcell,1) = CentroidMat(1,:);
            Results{k1}.CentroidY(startcell:endcell,1) = CentroidMat(2,:);
            Results{k1}.CytArea(startcell:endcell,1) = [cell2mat(CytArea)'; zeros(totcells_field-length(cell2mat(CytArea)),1)];
            
            % segment foci or subcellular structures
            if options.substructures.flag == 1 
                % Substructure channel
                i3 = options.substructures.channel;
                
                % Load Fluorsence Image without background
                FluorImage = imread(FileTif,'Index',i3);
                edgefluor = FluorImage == 0;
                edgeeliminatemask = stdfilt(edgefluor) == 0;
                FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*1)));
                Focisv_BK = FluorImage - imopen(FluorImage, strel('disk',options.substructures.bkgsize)); 
                SpatVar_Raw = stdfilt(Focisv_BK)./mean(FluorImage( FluorImage(:)>10  ));
                
                % Create mask based on Spatial Variation Threshold (1.0 default)
                FociSeg = SpatVar_Raw>options.substructures.SV_thresh;      
                FociSeg = FociSeg > 0 & edgeeliminatemask;
                FociSeg = imfill(FociSeg,'holes');

                % Dilate mask if needed
                Dilate_Nuc_Image = imdilate(lb_NucMask,offsetstrel('ball',options.substructures.nucdilation,0));
                lb_FociImage = uint16(FociSeg).*Dilate_Nuc_Image;
                
                % Fill in holes if needed
                lb_FociImage = imfill(lb_FociImage,'holes');
                mask = bwareaopen(lb_FociImage,10);
                lb_FociImage = lb_FociImage.*uint16(mask);
                FociSegImage = lb_FociImage > 0;

                % Remove Foci with odd shapes
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
    
                % Remove apoptotic cell by Area (removing rupture > 50% area)
                Fociarea = regionprops(lb_FociImage,"Area");
                Fociarea = [Fociarea.Area];
                Nucleararea = regionprops(Dilate_Nuc_Image,"Area");
                Nucleararea = [Nucleararea.Area];
                Fociarea(length(Fociarea)+1:length(Nucleararea)) = 0;
                Arearatio = Fociarea > (0.5 * Nucleararea);
                Arearatioindex = find(Arearatio);
                for i = Arearatioindex
                    lb_FociImage(lb_FociImage==i) = 0;
                end

                % Create Quality check files
                FociSegImage = lb_FociImage > 0;
                FociCheckImage = uint16(NucMask);
                FociCheckImage(FociSegImage>0)=2;
                DirectCheckIm = cat(3,uint16(FociCheckImage),uint16(Focisv_BK));
                DirectCheckIm = cat(3,DirectCheckIm,uint16(DAPI_img));
                substructuresegfile = [ilastikfolder 'checknucsubseg\' num2str(filename.realcoreinfo(k1).index) '_Seg_SubNucStruct_SV' num2str(options.substructures.SV_thresh) '.tif'];
                imwrite(uint16(DirectCheckIm),substructuresegfile)
    
                % Input Foci Information to Result Data
                Foci_stats = regionprops(lb_FociImage,FluorImage,'PixelValues');
                Foci_tot_stats = regionprops(Dilate_Nuc_Image,FluorImage,'PixelValues');
                Results{k1}.AreaSubstruct(1:max(lb_FociImage(:)),1) = cellfun(@length,{Foci_stats.PixelValues});
                Results{k1}.AreaSubstruct(1:max(Dilate_Nuc_Image(:)),2) = cellfun(@length,{Foci_tot_stats.PixelValues});
                Results{k1}.MeanInsideSubstruct(1:max(lb_FociImage(:)),1) = cellfun(@mean,{Foci_stats.PixelValues});
                Results{k1}.MeanInsideSubstruct(1:max(Dilate_Nuc_Image(:)),2) = cellfun(@mean,{Foci_tot_stats.PixelValues});
                Results{k1}.AreaSubstruct(isnan(Results{k1}.AreaSubstruct))=0;
                Results{k1}.MeanInsideSubstruct(isnan(Results{k1}.MeanInsideSubstruct))=0;

            end

            %%%% cycle around the images and extract measurements ---------
            for i3 = 1:max(filename.cycles)*4
                % Load each channel
                FluorImage = imread(FileTif,'Index',i3);
                FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));

                % Extract pixel values from fluo image using the segmented images    
                Nuclei_stats = regionprops(lb_NucMask,FluorImage_BK,'PixelValues');
                Cytopl_stats = regionprops(lb_CytMask,FluorImage_BK,'PixelValues');
                FullCell_stats = regionprops(lb_FullCellMask,FluorImage_BK,'PixelValues');
                Results{k1}.MeanNucSign(startcell:endcell,i3) = cellfun(@mean,{Nuclei_stats.PixelValues});
                Results{k1}.MeanCytSign(startcell:endcell,i3) = cellfun(@mean,{Cytopl_stats.PixelValues});
                Results{k1}.MeanFullCellSign(startcell:endcell,i3) = cellfun(@mean,{FullCell_stats.PixelValues});
    
                % All channel Foci segmentation
                Foci_stats = regionprops(lb_FociImage,FluorImage_BK,'PixelValues');
                Results{k1}.MeanFociSign(1:max(lb_FociImage(:)),i3) = cellfun(@mean,{Foci_stats.PixelValues});
                Results{k1}.MeanFociSign(isnan(Results{k1}.MeanFociSign))=0;
                Results{k1}.MeanFociSign( size(Results{k1}.Area,1),1) = 0;
            end
            % Mark as processed
            Results{k1}.CoreFlag = 1; 
            
            % Save data progress
            save([filename.resultfolder filename.resultfile],'Results', '-v7.3') %Saving matrix
        end
    end

    %Correction for NaN being deleted  
    for i = 1:length(Results)
        if ~isempty(Results{i})
            Results{i}.Corenum(1:size(Results{i}.Area,1),1) = i;
            if size(Results{i}.MeanFociSign,1) < size(Results{i}.Area,1) 
                Results{i}.MedianFociSign( size(Results{i}.Area,1),1) = uint16(0);
                Results{i}.MeanFociSign( size(Results{i}.Area,1),1) = 0;
            end
        end
    end
    % Save completed Result Data
    save([filename.resultfolder filename.resultfile],'Results', '-v7.3') %Saving matrix
end