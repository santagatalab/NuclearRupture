function Step4_CycIF_measurements_TMA(filename, options)  

    mkdir(filename.resultfolder); 
    mkdir([filename.microNucfolder filename.ilastiksegfol 'checknucsubseg\']);
    tic
    fullstackfolder = filename.fullstackfolder;
    ilastikfolder = [filename.microNucfolder filename.ilastiksegfol];
    ilastiksegfolder = [ilastikfolder 'seg\'];
    ilastikNucCytSegfolder = [ilastikfolder 'nuccytseg\'];
    mkdir(ilastikNucCytSegfolder);
    if exist([filename.resultfolder filename.resultfile],'file') == 2
        load([filename.resultfolder filename.resultfile])
    else
        Results{1} = [];
    end
    
    for k1 = 1:length(filename.realcoreinfo)
        totallength = length(filename.realcoreinfo);
        disp([num2str(k1/totallength *100) ' percent done']);
        FileTif = [fullstackfolder filename.realcoreinfo(k1).name];
        if ~isempty(FileTif)
            Results{k1} = [];
            Results{k1}.CoreFlag = [];
            Results{k1}.Area = [];
            Results{k1}.Solidity = [];
            Results{k1}.CentroidX = [];
            Results{k1}.CentroidY = [];
            Results{k1}.MeanNucSign = [];
            Results{k1}.MeanCytSign = [];
            %%%% load the FullStack and/to check it exists ----------------
    
            if exist(FileTif,'file')~=2
                disp(FileTif)
                disp('Error: FullStack not found')
                continue
            end
            disp(FileTif)
    
            %%%% load nuclear mask, create cytoplasmic mask and save them -
            %rootname = Results{k1}.Name;
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
    
    
            DAPI_img = uint16(imread(FileTif,'Index',1));
            NucCytSeg = cat(3,NucMask,DAPI_img,CytMask);
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
            %%%% cycle around the images and extract measurements ---------
            for i3 = 1:max(filename.cycles)*4
                % load image
                FluorImage = imread(FileTif,'Index',i3);
                FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
                % extract pixel values from fluo image using the segmented images    
                Nuclei_stats = regionprops(lb_NucMask,FluorImage_BK,'PixelValues');
                Cytopl_stats = regionprops(lb_CytMask,FluorImage_BK,'PixelValues');
                FullCell_stats = regionprops(lb_FullCellMask,FluorImage_BK,'PixelValues');
    
                Results{k1}.MeanNucSign(startcell:endcell,i3) = cellfun(@mean,{Nuclei_stats.PixelValues});
                Results{k1}.MeanCytSign(startcell:endcell,i3) = cellfun(@mean,{Cytopl_stats.PixelValues});
                Results{k1}.MeanFullCellSign(startcell:endcell,i3) = cellfun(@mean,{FullCell_stats.PixelValues});
            end
            Results{k1}.CoreFlag = 1; 
            save([filename.resultfolder filename.resultfilemicro],'Results', '-v7.3') %Saving matrix
        end
    end


    
end
