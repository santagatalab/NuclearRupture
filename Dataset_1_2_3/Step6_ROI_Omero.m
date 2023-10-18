function Step6_ROI_Omero(filename,options)
    % Coordinates correction and ROI mask creation
    load([filename.resultfolder filename.resultfile])
    save([filename.resultfolder 'Step5backup_' filename.resultfile],'Results', '-v7.3')
    for coren = 1:length(Results)
        if isempty(Results{coren})
            continue
        end
        % Extract fieldname rows and columns
        cName = filename.realcoreinfo(coren).name;
        [~,cName_seg] = strtok(cName,'_');
        [~,cName_seg2] = strtok(cName_seg,'_');
        [yn_string,cName_seg3] = strtok(cName_seg2(2:end),'_');
        [xn_string,~] = strtok(cName_seg3(2:end),'.tif');
        xn = str2num(xn_string);
        yn = str2num(yn_string);
        for i = 1:length(Results{coren}.Area)
             Results{coren}.FullCentroidX(i,1) = Results{coren}.CentroidX(i) + 5000*(xn-1);
             Results{coren}.FullCentroidY(i,1)= Results{coren}.CentroidY(i) + 5000*(yn-1);
        end
    end
    disp('Coordinates correction done')

    % Read ROI file from Omero
    T = readtable(filename.roifile);
    % Size of the image
    filenamesv = [filename.datafolder filename.folders.coordinates filename.tissues{1}];
    load(filenamesv)
    bottomrightcoords = Coordinates.Field{end,end};
    maxy = bottomrightcoords(5);
    maxx = bottomrightcoords(3);
    % Create ROI mask
    ROImask = zeros(maxx,maxy);
    k_changing = 0;
    regionroi =[];
    for k = 1:length(T.Points)
        if k > 1 && strcmp(T.Points{k},T.Points{k-1})
            continue
        end
        currentTpoints = T.Points{k};
        croilist_x = [];
        croilist_y = [];
        while ~isempty(currentTpoints)
            [cxvalue, currentTpoints] = strtok(currentTpoints,',');
            currentTpoints = currentTpoints(2:end);
            [cyvalue, currentTpoints] = strtok(currentTpoints,' ');
            currentTpoints = currentTpoints(2:end);
            croilist_x = [croilist_x str2num(cxvalue)];
            croilist_y = [croilist_y str2num(cyvalue)];
        end
        if isempty(T.text{k}) %general ROI
            k_changing = k_changing+1; 
            xroilist{k_changing} = croilist_x;
            yroilist{k_changing} = croilist_y;
        else %specific ROI
            currentregionname = T.text{k};
            % Remove illegal characters
            currentregionname(currentregionname == ' ') = [];
            currentregionname(currentregionname == '(') = '_';
            currentregionname(currentregionname == ')') = [];
            currentregionname(currentregionname == '.') = '_';
            currentregionname(currentregionname == '-') = '_';
            % Assign ROI points information to structure regionroi
            if ~exist(['regionroi.' currentregionname ''],"var")
                regionroi.(currentregionname) = [];
                regionroi.(currentregionname).x{1} = croilist_x;
                regionroi.(currentregionname).y{1} = croilist_y;
                regionroi.(currentregionname).kc = 1;
            else
                kc = regionroi.(currentregionname).kc;
                regionroi.(currentregionname).x{kc+1} = croilist_x;
                regionroi.(currentregionname).y{kc+1} = croilist_y;
                regionroi.(currentregionname).kc = kc + 1;
            end
        end
    end

    %
    Regionnames = fieldnames(regionroi);
    for coren = 1:length(Results)
        if isempty(Results{coren})
            continue
        end
        Results{coren}.roimask = [];
        for celln = 1:length(Results{coren}.Area)
            cx = Results{coren}.FullCentroidX(celln);
            cy = Results{coren}.FullCentroidY(celln);
            Results{coren}.regionassign{celln} = [];
            Results{coren}.roimask(celln) = false;
            for k = 1:length(xroilist)
                if ~Results{coren}.roimask(celln) && inpolygon(cx,cy,xroilist{k},yroilist{k})
                    Results{coren}.roimask(celln) = true;
                end
            end
            for regionroi_nameindex = 1:length(Regionnames)
                currentregionname = Regionnames{regionroi_nameindex};
                for k = 1:regionroi.(currentregionname).kc
                    croilist_x = regionroi.(currentregionname).x{k};
                    croilist_y = regionroi.(currentregionname).y{k};
                    if ~inpolygon(cx,cy,croilist_x,croilist_y)
                        continue
                    end
                    if isempty(Results{coren}.regionassign{celln}) 
                        Results{coren}.regionassign{celln} = {currentregionname};
                    else
                        Results{coren}.regionassign{celln} = [Results{coren}.regionassign{celln} {currentregionname}];
                    end
                end                
            end
        end
        Results{coren}.roimask = Results{coren}.roimask';
        Results{coren}.regionassign = Results{coren}.regionassign';
    end
    %Remove everything that is not in the ROI
    Fname =  fieldnames(Results{1});
    newResults = [];
    for coren = 1:length(Results)
        if isempty(Results{coren})
            continue
        end
        cresult = Results{coren}.Area;
        [cmaxx_area,cmaxy_area] = size(cresult);
        for Fnamen = 1:length(Fname)
            cFname = Fname{Fnamen};
            if strcmp(cFname,'roimask')
                continue
            end
            if isempty(Results{coren}.Area)
                continue
            end
            cresult = Results{coren}.(cFname);
            [cmaxx,cmaxy] = size(cresult);
            if cmaxx ~= cmaxx_area
                newResults{coren}.(cFname) = cresult;
            else
                newResults{coren}.(cFname) = cresult((Results{coren}.roimask==1),:);
            end
        end
    end
    Results = newResults;

    % Narrow down to only one assignment for each cell
     for coren = 1:length(Results)
        if isempty(Results{coren})
            continue
        end
        Results{coren}.tumor = [];
        Results{coren}.type = [];
        for celln = 1:length(Results{coren}.Area)
            Results{coren}.tumor(celln) = false;
           if isempty(Results{coren}.regionassign{celln})
               Results{coren}.type{celln} = 'NoType';
               continue
           end
           currentRegionAssign = Results{coren}.regionassign{celln};
           tumor_bulkindex = [];
           % If assignments contain tumor_bulk, classify as tumor
           for i = 1:length(currentRegionAssign)
               if contains(currentRegionAssign{i},'Tumor_Bulk')
                    Results{coren}.tumor(celln) = true;
                    tumor_bulkindex = i;
                    Results{coren}.type{celln} = 'Tumor';
               end
           end
           if ~isempty(tumor_bulkindex)
               currentRegionAssign(tumor_bulkindex) = [];
           end
           if length(currentRegionAssign) == 1
               Results{coren}.type{celln} = currentRegionAssign{1};
           elseif length(currentRegionAssign) > 1
               Results{coren}.type{celln} = 'MultiType';
           end
        end
        Results{coren}.tumor = Results{coren}.tumor';
        Results{coren}.type = Results{coren}.type';
     end 
    save([filename.resultfolder filename.resultfile],'Results', '-v7.3') %Saving matrix

    %%Micronucleus
    % Check if micronucleus result file exists
    if exist([filename.resultfolder filename.resultfilemicro],'file')
        % Coordinates correction and ROI mask creation
        clear Results
        load([filename.resultfolder filename.resultfilemicro])
        save([filename.resultfolder 'Step5backup_' filename.resultfilemicro],'Results', '-v7.3')
        for coren = 1:length(Results)
            if isempty(Results{coren})
                continue
            end
            cName = filename.realcoreinfo(coren).name;
            [~,cName_seg] = strtok(cName,'_');
            [~,cName_seg2] = strtok(cName_seg,'_');
            [yn_string,cName_seg3] = strtok(cName_seg2(2:end),'_');
            [xn_string,~] = strtok(cName_seg3(2:end),'.tif');
            xn = str2num(xn_string);
            yn = str2num(yn_string);
            for i = 1:length(Results{coren}.Area)
                 Results{coren}.FullCentroidX(i,1) = Results{coren}.CentroidX(i) + 5000*(xn-1);
                 Results{coren}.FullCentroidY(i,1)= Results{coren}.CentroidY(i) + 5000*(yn-1);
            end
        end
        disp('Coordinates correction done')
    
        % Read ROI file
        T = readtable(filename.roifile);
        % Size of the image
        filenamesv = [filename.datafolder filename.folders.coordinates filename.tissues{1}];
        load(filenamesv)
        bottomrightcoords = Coordinates.Field{end,end};
        maxy = bottomrightcoords(5);
        maxx = bottomrightcoords(3);
        % Create ROI mask
        ROImask = zeros(maxx,maxy);
        k_changing = 0;
        regionroi =[];
        for k = 1:length(T.Points)
            if k > 1 && strcmp(T.Points{k},T.Points{k-1})
                continue
            end
            currentTpoints = T.Points{k};
            croilist_x = [];
            croilist_y = [];
            while ~isempty(currentTpoints)
                [cxvalue, currentTpoints] = strtok(currentTpoints,',');
                currentTpoints = currentTpoints(2:end);
                [cyvalue, currentTpoints] = strtok(currentTpoints,' ');
                currentTpoints = currentTpoints(2:end);
                croilist_x = [croilist_x str2num(cxvalue)];
                croilist_y = [croilist_y str2num(cyvalue)];
            end
            if isempty(T.text{k}) %general ROI
                k_changing = k_changing+1; 
                xroilist{k_changing} = croilist_x;
                yroilist{k_changing} = croilist_y;
            else %specific ROI
                currentregionname = T.text{k};
                currentregionname(currentregionname == ' ') = [];
                currentregionname(currentregionname == '(') = '_';
                currentregionname(currentregionname == ')') = [];
                currentregionname(currentregionname == '.') = '_';
                currentregionname(currentregionname == '-') = '_';
                if ~exist(['regionroi.' currentregionname ''],"var")
                    regionroi.(currentregionname) = [];
                    regionroi.(currentregionname).x{1} = croilist_x;
                    regionroi.(currentregionname).y{1} = croilist_y;
                    regionroi.(currentregionname).kc = 1;
                else
                    kc = regionroi.(currentregionname).kc;
                    regionroi.(currentregionname).x{kc+1} = croilist_x;
                    regionroi.(currentregionname).y{kc+1} = croilist_y;
                    regionroi.(currentregionname).kc = kc + 1;
                end
            end
        end
        Regionnames = fieldnames(regionroi);
        for coren = 1:length(Results)
            if isempty(Results{coren})
                continue
            end
            Results{coren}.roimask = [];
            for celln = 1:length(Results{coren}.Area)
                cx = Results{coren}.FullCentroidX(celln);
                cy = Results{coren}.FullCentroidY(celln);
                Results{coren}.regionassign{celln} = [];
                Results{coren}.roimask(celln) = false;
                for k = 1:length(xroilist)
                    if ~Results{coren}.roimask(celln) && inpolygon(cx,cy,xroilist{k},yroilist{k})
                        Results{coren}.roimask(celln) = true;
                    end
                end
                for regionroi_nameindex = 1:length(Regionnames)
                    currentregionname = Regionnames{regionroi_nameindex};
                    for k = 1:regionroi.(currentregionname).kc
                        croilist_x = regionroi.(currentregionname).x{k};
                        croilist_y = regionroi.(currentregionname).y{k};
                        if ~inpolygon(cx,cy,croilist_x,croilist_y)
                            continue
                        end
                        if isempty(Results{coren}.regionassign{celln}) 
                            Results{coren}.regionassign{celln} = {currentregionname};
                        else
                            Results{coren}.regionassign{celln} = [Results{coren}.regionassign{celln} {currentregionname}];
                        end
                    end                
                end
            end
            Results{coren}.roimask = Results{coren}.roimask';
        end
        %Remove everything that is not in the ROI
        Fname =  fieldnames(Results{2});
        newResults = [];
        for coren = 1:length(Results)
            if isempty(Results{coren})
                continue
            end
            cresult = Results{coren}.Area;
            [cmaxx_area,cmaxy_area] = size(cresult);
            for Fnamen = 1:length(Fname)
                cFname = Fname{Fnamen};
                if strcmp(cFname,'roimask')
                    continue
                end
                if isempty(Results{coren}.Area)
                    continue
                end
                cresult = Results{coren}.(cFname);
                [cmaxx,cmaxy] = size(cresult);
                if cmaxx ~= cmaxx_area
                    newResults{coren}.(cFname) = cresult;
                else
                    newResults{coren}.(cFname) = cresult((Results{coren}.roimask==1),:);
                end
            end
        end
        Results = newResults;
        %%%narrow down to only one assignment
        for coren = 1:length(Results)
            if isempty(Results{coren})
                continue
            end
            Results{coren}.tumor = [];
            Results{coren}.type = [];
            for celln = 1:length(Results{coren}.Area)
                Results{coren}.tumor(celln) = false;
               if isempty(Results{coren}.regionassign{celln})
                   Results{coren}.type{celln} = 'NoType';
                   continue
               end
               currentRegionAssign = Results{coren}.regionassign{celln};
               tumor_bulkindex = [];
               for i = 1:length(currentRegionAssign)
                   if contains(currentRegionAssign{i},'Tumor_Bulk')
                        Results{coren}.tumor(celln) = true;
                        tumor_bulkindex = i;
                        Results{coren}.type{celln} = 'Tumor';
                   end
               end
               if ~isempty(tumor_bulkindex)
                   currentRegionAssign(tumor_bulkindex) = [];
               end
               if length(currentRegionAssign) == 1
                   Results{coren}.type{celln} = currentRegionAssign{1};
               elseif length(currentRegionAssign) > 1
                   Results{coren}.type{celln} = 'MultiType';
               end
            end
            Results{coren}.tumor = Results{coren}.tumor';
            Results{coren}.type = Results{coren}.type';
        end 
        save([filename.resultfolder filename.resultfilemicro],'Results', '-v7.3') %Saving matrix
    end
end

