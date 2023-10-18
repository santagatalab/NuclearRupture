function Analysis_TumorBoundary(filename,options)
    % Read ROI file
    T = readtable(filename.roifile);
    % Size of the image
    filenamesv = [filename.datafolder filename.folders.coordinates filename.tissue];
    try
        mkdir(filename.graphsfolder);
    end
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

    TumorBoundary = [regionroi.Tumor_Bulk.x{:};regionroi.Tumor_Bulk.y{:} ];
    %TumorBoundary = [regionroi.Tumor_Peri_necrotic.x{:};regionroi.Tumor_Peri_necrotic.y{:} ];
    %% Load aggregate data
      load([filename.resultfolder 'Aggregate_' filename.resultfile])

    Channelcutoffcell = readcell([filename.datafolder 'Channelscutoff_BC.xlsx']);
    [maxchannel,~] = size(Channelcutoffcell);
    AllChannel = Channelcutoffcell(:,1);
    AllChanneln = 1:maxchannel;
    Channelcutoff = cell2mat(Channelcutoffcell(:,2));
    amps = cell2mat(Channelcutoffcell(:,3));
    Channelchar = Channelcutoffcell(:,4);
    ChannelL = Channelcutoffcell(:,5);
    ChannelC = Channelcutoffcell(:,6);
    %% Add tumor boundary
    AggrResults_primary.TumorBulkDist = [];
    for i = 1:length(AggrResults_primary.Area)
        cx = AggrResults_primary.FullCentroidX(i);
        cy = AggrResults_primary.FullCentroidY(i);
        [d,~,~] = p_poly_dist(cx, cy, TumorBoundary(1,:), TumorBoundary(2,:));
        AggrResults_primary.TumorBulkDist(i,1) = d;
    end

    %% Graph
    a = figure();
    clist = [48 60 50 56 18]; % 48 60 50
    lineclourlist = {"#EDB120",'r','b','cyan','green'};
    
    legendnamelist = cell(1,length(clist)*2);
    for c_index = 1:length(clist)
        c = clist(c_index);
        cchannelname = AllChannel{c};
        pixelsize = 0.325;

        x_unique = [];
        y_unique = [];
        y_avg = [];
        y_top75 = [];
        y_bot25 = [];
        y_std = [];
        x = AggrResults_primary.TumorBulkDist * pixelsize;
        y = log10(AggrResults_primary.MeanFullCellSign(:,c));
        [x,sortedindex] = sort(x,'ascend');
        y = y(sortedindex);
        x = round(x,-2);
        [x_unique,~,ib] = unique(x);
%         [ib,edges] = discretize(x,500);
%         x_unique = 
        for i = 1:length(x_unique)
            y_unique{i} = y(ib==i);
            validindex=y~=-Inf;
            y_avg(i) = nanmean(y(ib==i & validindex));
            y_top75(i) = prctile(y_unique{i},75);
            y_bot25(i) = prctile(y_unique{i},25);
            y_std(i) = std(y_unique{i});
        end
        % Remove small samples
        delmask = x_unique >0;
        for i = 1:length(x_unique)
            delmask(i) = length(y_unique{i}) < 300;
        end
        x_unique(delmask) = [];
        y_unique(delmask) = [];
        y_avg(delmask) = [];
        y_top75(delmask) = [];
        y_bot25(delmask) = [];
        y_std(delmask) = [];

        %scatter(x,y,1)
        %s = scatter(x_unique,y_avg,5);
        x_line = x_unique;
        y_line = smoothdata(y_avg,"SmoothingFactor",0.1);
        s = plot(x_line,y_line,'LineWidth',2,'Color',lineclourlist{c_index});
        hold on
        %shade data rage
         t  = x_unique';
         y1 = y_bot25;
         y2 = y_top75;
    %     y1 = y1 - y_std;
    %     y2 = y2 + y_std;
         y1 = smoothdata(y1,"SmoothingFactor",0.1);
         y2 = smoothdata(y2,"SmoothingFactor",0.1);
         f = fill([t fliplr(t)],[y1 fliplr(y2)] ,'r','EdgeColor','none');
         f.FaceColor = lineclourlist{c_index};
         alpha(f,0.1)
         legendnamelist((c_index-1)*2+1) = {cchannelname};
         legendnamelist((c_index-1)*2+2) = {''};
    end

    maxy = 5;
    xlim([-3000 2000])
    %xlim([min(x_unique)*1.1 max(x_unique)*1.1])
    %ylim([0 maxy])
    xline(0,'--');
    text(0,maxy*0.98,' Outside');
    text(0,maxy*0.98,'Inside ','HorizontalAlignment', 'right');
    legendnamelist(end+1) = {''};
    legend(legendnamelist)
    xlabel('Distance from Tumor Boundary (µm)')
    ylabel('Log10(Mean Nuclear Intensity)')
    title(['Spacial Distribution of Lineage Markers in Autopsy ' filename.tissues{1}])
    exportgraphics(a,[filename.graphsfolder 'Spacial Distribution of Lineage Markers in Autopsy ' filename.tissues{1},'.pdf'],'Resolution',600,'ContentType','vector')
    %% Graph
    b = figure();
    clist = [20 46 54]; % 48 60 50
    lineclourlist = {"cyan",'magenta','r'};
    
    legendnamelist = cell(1,length(clist)*2);
    for c_index = 1:length(clist)
        c = clist(c_index);
        cchannelname = AllChannel{c};
        pixelsize = 0.325;

        x_unique = [];
        y_unique = [];
        y_avg = [];
        y_top75 = [];
        y_bot25 = [];
        y_std = [];
        x = AggrResults_primary.TumorBulkDist * pixelsize;
        y = log10(AggrResults_primary.MeanFullCellSign(:,c));
        [x,sortedindex] = sort(x,'ascend');
        y = y(sortedindex);
        x = round(x,-2);
        [x_unique,~,ib] = unique(x);
%         [ib,edges] = discretize(x,500);
%         x_unique = 
        for i = 1:length(x_unique)
            y_unique{i} = y(ib==i);
            validindex=y~=-Inf;
            y_avg(i) = nanmean(y(ib==i & validindex));
            y_top75(i) = prctile(y_unique{i},75);
            y_bot25(i) = prctile(y_unique{i},25);
            y_std(i) = std(y_unique{i});
        end
        % Remove small samples
        delmask = x_unique >0;
        for i = 1:length(x_unique)
            delmask(i) = length(y_unique{i}) < 300;
        end
        x_unique(delmask) = [];
        y_unique(delmask) = [];
        y_avg(delmask) = [];
        y_top75(delmask) = [];
        y_bot25(delmask) = [];
        y_std(delmask) = [];

        %scatter(x,y,1)
        %s = scatter(x_unique,y_avg,5);
        y_line = smoothdata(y_avg,"SmoothingFactor",0.1);
        x_line = x_unique;
        s = plot(x_line,y_line,'LineWidth',2,'Color',lineclourlist{c_index});
        hold on
        %shade data rage
         t  = x_unique';
         y1 = y_bot25;
         y2 = y_top75;
    %     y1 = y1 - y_std;
    %     y2 = y2 + y_std;
         y1 = smoothdata(y1,"SmoothingFactor",0.1);
         y2 = smoothdata(y2,"SmoothingFactor",0.1);
         f = fill([t fliplr(t)],[y1 fliplr(y2)] ,'r','EdgeColor','none');
         f.FaceColor = lineclourlist{c_index};
         alpha(f,0.1)
         legendnamelist((c_index-1)*2+1) = {cchannelname};
         legendnamelist((c_index-1)*2+2) = {''};
    end

    maxy = 3;
    xlim([-3000 2000])
    %xlim([min(x_unique)*1.1 max(x_unique)*1.1])
    %ylim([0 maxy])
    xline(0,'--');
    text(0,maxy*0.98,' Outside');
    text(0,maxy*0.98,'Inside ','HorizontalAlignment', 'right');
    legendnamelist(end+1) = {''};
    legend(legendnamelist)
    xlabel('Distance from Tumor Boundary (µm)')
    ylabel('Log10(Mean Nuclear Intensity)')
    title(['Spacial Distribution of BAF and Lamins Markers in Autopsy ' filename.tissues{1}])
    exportgraphics(b,[filename.graphsfolder 'Spacial Distribution of BAF and Lamins Markers in Autopsy ' filename.tissues{1},'.pdf'],'Resolution',600,'ContentType','vector')

    %% Test
    figure
    pgon = polyshape(TumorBoundary(1,:),TumorBoundary(2,:));
    plot(pgon)
end
