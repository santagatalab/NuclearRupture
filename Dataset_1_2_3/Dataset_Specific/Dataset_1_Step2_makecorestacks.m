function Dataset_1_Step2_makecorestacks(filename,options)

realcorefolder =  [filename.analfolder 'RealCores' filesep];
fullstackfolder = [filename.analfolder 'FullStacks' filesep];
dapistackfolder = [filename.analfolder 'DapiStacks' filesep];
mkdir(fullstackfolder)
mkdir(dapistackfolder)

half = options.tilesize/2;

% cycle through the coresng 
% for i1 = 1:length(filename.prefix1)
%     for i2 = 1:length(filename.prefix2)

for i1 = 1:length(filename.realcoreinfo)
    if ~isempty(filename.realcoreinfo{i1})
        
        flag = [];
        
        final_stack = [fullstackfolder 'Core'      filename.realcoreinfo{i1}.name filename.suffix ];
        DAPI_stack  = [dapistackfolder 'DAPI_Core' filename.realcoreinfo{i1}.name filename.suffix ];
        DAPIStack = [];
        RunningStack = [];
        
        img_name = [realcorefolder 'RealCore_' filename.realcoreinfo{i1}.name filename.suffix];
        try
            DAPICycle0 = imread(img_name);
        catch
            disp(['MakeCoreStacks Error: RealCore_' filename.realcoreinfo{i1}.name ' not found'])
            continue
        end
        
        if min(size(DAPICycle0)) < half
            disp(['MakeCoreStacks Error: RealCore ' filename.realcoreinfo{i1}.name ' is too small'])
            continue
        end
        
        if exist(final_stack,'file') == 2
            disp(['FullStack already made for ' filename.realcoreinfo{i1}.name])
            continue
        end

%         try
%             for i = 1:options.num_old_DAPI
%                 DAPIStack(:,:,i) = imread(DAPI_stack,'Index',i);
%                 RunningStack(:,:,i) = imread(final_stack,'Index',i);
%             end
%         catch
            RunningStack = [];
            DAPIStack = [];
%         end
        
        
        
        
        
        for i3 = 1:length(filename.cycles) 
            disp(['Started Core ' filename.realcoreinfo{i1}.name ' cycle ' num2str(filename.cycles(i3))])
            clear mont_row mont_col tile_row tile_col flag
            % initialize new DAPI images with all zeros
            DAPI{i3} = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');
            
            % load the DAPI images from the next cycle
            for j1 = 1:length(filename.realcoreinfo{i1}.tiles)
                
                tiles_img{j1} = imread([filename.datafolder filename.cycleprefix num2str(filename.cycles(i3)) filesep ...
                                        filename.realcoreinfo{i1}.name '(fld ' num2str(filename.realcoreinfo{i1}.tiles(j1),filename.dim2) ...
                                        filename.wavelengths{1} filename.suffix]);
                
                                                                       

                    for j2 = 1:4
                        tiles_part{j1,j2} = tiles_img{j1}((1:half)+half*floor(j2/3),(1:half)+half*(1-mod(j2,2)));
                        
                        c = normxcorr2(tiles_part{j1,j2},DAPICycle0);

                    
                        [ypeak, xpeak] = find(c==max(c(:)));
                        tileshift_y(j1,j2) = ypeak-size(tiles_part{j1,j2},1);
                        tileshift_x(j1,j2) = xpeak-size(tiles_part{j1,j2},2);
                        edgestrength(j1,j2) = max(c(:));
                        if max(c(:)) > options.corrthresh
                            flag{j1,j2} = 1;
                            mont_row{j1,j2} = 1+max(0,tileshift_y(j1,j2)):min(tileshift_y(j1,j2)+size(tiles_part{j1,j2},1),size(DAPI{i3},1));
                            mont_col{j1,j2} = 1+max(0,tileshift_x(j1,j2)):min(tileshift_x(j1,j2)+size(tiles_part{j1,j2},2),size(DAPI{i3},2));
                            tile_row{j1,j2} = (1-min(0,tileshift_y(j1,j2))):(size(tiles_part{j1,j2},1)-max(0,tileshift_y(j1,j2)+size(tiles_part{j1,j2},1)-size(DAPI{i3},1)));
                            tile_col{j1,j2} = (1-min(0,tileshift_x(j1,j2))):(size(tiles_part{j1,j2},2)-max(0,tileshift_x(j1,j2)+size(tiles_part{j1,j2},2)-size(DAPI{i3},2)));
                
                    
                            DAPI{i3}(mont_row{j1,j2},mont_col{j1,j2}) = tiles_part{j1,j2}(tile_row{j1,j2},tile_col{j1,j2}); 
%                             figure
%                             imshow(DAPI{i3},[])
%                             disp(num2str(max(c(:))))
                            
                        else
                            flag{j1,j2} = 0;
%                             disp(['MakeCoreStacks Error: Rejected a tile with max(c(:)) = ' num2str(max(c(:))) ])
                        end
                    end
                    
            end 

        DAPIStack = cat(3,DAPIStack,DAPI{i3});
        RunningStack = cat(3,RunningStack,DAPI{i3});
        
        
        for wv1 = 2:length(filename.wavelengths)
            temp_mont = zeros(size(DAPICycle0,1),size(DAPICycle0,2),'uint16');
            

            for j1 = 1:length(filename.realcoreinfo{i1}.tiles)

                tiles_img{j1} = imread([filename.datafolder filename.cycleprefix num2str(filename.cycles(i3)) filesep ...
                                        filename.realcoreinfo{i1}.name '(fld ' num2str(filename.realcoreinfo{i1}.tiles(j1),filename.dim2) ...
                                        filename.wavelengths{wv1} filename.suffix]);

                                        
                for j2 = 1:4
                    if flag{j1,j2} == 1
                        tiles_part{j1,j2} = tiles_img{j1}((1:half)+half*floor(j2/3),(1:half)+half*(1-mod(j2,2)));
                        temp_mont(mont_row{j1,j2},mont_col{j1,j2}) = tiles_part{j1,j2}(tile_row{j1,j2},tile_col{j1,j2}); 

                    end
                end
            end
            
            RunningStack = cat(3,RunningStack,temp_mont);
            clear  temp_mont_cropped temp_mont tiles_img
        end

            
        end
               
    optionssave.append = options.appendflag;
    saveastiff(RunningStack,final_stack,optionssave)
    saveastiff(DAPIStack,DAPI_stack,optionssave)
    clear RunningStack DAPIStack DAPI_img DAPICycle0
        
    end
end


