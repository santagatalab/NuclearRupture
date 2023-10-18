function Dataset_1_Step1_RealCoreCreation(filename,options)

slice = options.tilesize/options.scale;
maxtile = options.rows*options.cols;

% not to input - general parameters
uporlx = 1:slice;
dnorrx = (options.tilesize-slice+1):options.tilesize;
edgematrix = [1 1 0; 2 2 3; 3 3 4; 4 1 0];

count = 0;

for i1 = 1:length(filename.prefix1)
    for i2 = 1:length(filename.prefix2)
        
        realcorefile = [filename.analfolder 'RealCores' filesep 'RealCore_' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),options.dim) filename.suffix];
        disp(['Started ' 'RealCore_' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),options.dim)])
        if exist(realcorefile,'file') == 2
            count = count + 1;
            continue  
        end
        
        
        try
        % load the 4 files
        for j1 = 1:maxtile
            tilename = [filename.datafolder filename.realcorecycleprefix filesep filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),options.dim) '(fld ' num2str(options.fileordchange(j1)) filename.wavelengths{1} filename.suffix];
            tiles_img{j1} = imread(tilename);
            test(j1) = prctile(tiles_img{j1}(:),99);
        end
        catch
            disp(['Core' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),options.dim) ' is not found'])
            disp(tilename)
            continue
        end
        
        if max(test) < options.background
            disp(['Core' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),options.dim) ' might be empty'])
            continue
        else
        end
                
        % calcuate the cross corrs for the 4 edges
        
        try
            c{1} = normxcorr2(locallapfilt(tiles_img{1}(:,dnorrx),options.sigma,options.alpha),locallapfilt(tiles_img{2}(:,uporlx),options.sigma,options.alpha));
        catch
            c{1} = normxcorr2(tiles_img{1}(:,dnorrx),tiles_img{2}(:,uporlx));
        end
        try
            c{2} = normxcorr2(locallapfilt(tiles_img{2}(dnorrx,:),options.sigma,options.alpha),locallapfilt(tiles_img{3}(uporlx,:),options.sigma,options.alpha)); 
        catch
            c{2} = normxcorr2(tiles_img{2}(dnorrx,:),tiles_img{3}(uporlx,:)); 
        end
        try
            c{3} = normxcorr2(locallapfilt(tiles_img{3}(:,uporlx),options.sigma,options.alpha),locallapfilt(tiles_img{4}(:,dnorrx),options.sigma,options.alpha));
        catch
            c{3} = normxcorr2(tiles_img{3}(:,uporlx),tiles_img{4}(:,dnorrx));
        end
        try           
            c{4} = normxcorr2(locallapfilt(tiles_img{1}(dnorrx,:),options.sigma,options.alpha),locallapfilt(tiles_img{4}(uporlx,:),options.sigma,options.alpha));
        catch
            c{4} = normxcorr2(tiles_img{1}(dnorrx,:),tiles_img{4}(uporlx,:));
        end
        
        % calculate the max and the offset for each edge
        for j1 = 1:4
            edgestrength(j1) = max(c{j1}(:));
            [rowpeak, colpeak] = find(c{j1}==max(c{j1}(:)));
            if length(rowpeak)==1
                offset{j1}(1) = rowpeak-(size(c{j1},1)+1)/2;
            else
                offset{j1}(1) = (size(c{j1},1)+1)/2;
            end
            if length(colpeak)==1
                offset{j1}(2) = colpeak-(size(c{j1},2)+1)/2;
            else
                offset{j1}(1) = (size(c{j1},1)+2)/2;
            end
        end
    
        % determine path
        path(1,:) = [1 0 0];
        if edgestrength(1) > edgestrength(4)
            % use edge 1 to place tile 2 wrt 1
            path(2,:) = [2 1 1];
            % and move on
            if edgestrength(2) > edgestrength(4)
                path(3,:) = [3 2 1];
                if edgestrength(3) > edgestrength(4)
                    path(4,:) = [4 3 1];
                else
                    path(4,:) = [4 4 1];
                end 
            else
                path(3,:) = [4 4 1];
                if edgestrength(2) > edgestrength(3)
                    path(4,:) = [3 2 1];
                else
                    path(4,:) = [3 3 2];
                end 
            end
        else
            % use edge 4 to place tile 4 wrt 1
            path(2,:) = [4 4 1];
            % and move on
            if edgestrength(1) > edgestrength(3)
                path(3,:) = [2 1 1];
                if edgestrength(2) > edgestrength(3)
                    path(4,:) = [3 2 1];
                else
                    path(4,:) = [3 3 2];
                end 
            else
                path(3,:) = [3 3 2];
                if edgestrength(1) > edgestrength(2)
                    path(4,:) = [2 1 1];
                else
                    path(4,:) = [2 2 2];
                end 
            end
        end
        
        % define un-registered positions assuming 0 shift

        pos{1} = [0 0];
        pos{2} = [];
        pos{3} = [];
        pos{4} = [];
        
        % update the position based on the shifts
        for j1 = 2:maxtile
            edgeinuse = path(j1,2);
            % flip the offset in case the stitching direction is not in the
            % right one
            if path(j1,3) == 1
                flip = -1;
            else
                flip = 1;
            end
            
            pos{path(j1,1)} = pos{edgematrix(path(j1,2),path(j1,3)+1)} + offset{edgeinuse}*flip;          
        end
        
        pos{1} = [options.buffer options.buffer] + pos{1} + [0 0];
        pos{2} = [options.buffer options.buffer] + pos{2} + [0 options.tilesize-slice];
        pos{3} = [options.buffer options.buffer] + pos{3} + [options.tilesize-slice options.tilesize-slice];
        pos{4} = [options.buffer options.buffer] + pos{4} + [options.tilesize-slice 0];

        
        
        % build the montage
        
        DAPI_mont = zeros((options.tilesize+options.buffer)*2,(options.tilesize+options.buffer)*2);
        for j1 = 1:maxtile
            j2 = path(5-j1,1);
            if pos{j2}(1) < 0 || pos{j2}(2) < 0
                continue
            else
                DAPI_mont(pos{j2}(1)+1:pos{j2}(1)+options.tilesize,pos{j2}(2)+1:pos{j2}(2)+options.tilesize) = tiles_img{j2};
            end
        end
      
        % shave off the edges
        
        opened_img = imclose(DAPI_mont,strel('disk',200));
        core_outline = opened_img > options.background;
        
        core_outline = bwareaopen(core_outline,3.15*(options.coreradius^2));
        stats = regionprops(core_outline,'Centroid','Area','BoundingBox');
        
        area = cat(1,stats.Area);
        
        if length(area)==0
            disp(['Could not detect a core in Core' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),options.dim)])
            continue
        else


            [val,ind] = max(area);
            if length(area)>1
                disp(['Core' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),options.dim) 'might not have worked'])
            end

            col_start = round(stats(ind).BoundingBox(1));
            col_end   = col_start+round(stats(ind).BoundingBox(3))-1;
            row_start = round(stats(ind).BoundingBox(2));
            row_end   = row_start+round(stats(ind).BoundingBox(4))-1;

            row_core = row_start:row_end;
            col_core = col_start:col_end;
            realcore = DAPI_mont(row_core,col_core);

            options.message = 'false';
            saveastiff(uint16(realcore),[filename.analfolder filesep 'RealCores\RealCore_' filename.prefix1{i1} ' - '  num2str(filename.prefix2(i2),options.dim) filename.suffix],options)     
        end
        
    end
end
