clear
basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF';
savedir = [basedir '\Analysis\gfp'];
segdir = [basedir '\Analysis\seg'];

mkdir(segdir)

% index = 10;
% filename = [savedir '\' num2str(index) '_gfp.tiff'];
% tiff_info = imfinfo(filename);
% gfp_stack = imread(filename,1);
% %concatenate each successive tiff to tiff_stack
% for ii = 2 : size(tiff_info, 1)
%     temp_tiff = imread(filename, ii);
%     gfp_stack = cat(3 , gfp_stack, temp_tiff);
% end
% clear ii 

for index = 1:10
filename = [savedir '\' num2str(index) '_gfp.tiff'];
tiff_info = imfinfo(filename);
gfp_stack = imread(filename,1);
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(filename, ii);
    gfp_stack = cat(3 , gfp_stack, temp_tiff);
end
clear ii 
gpf_cell{index} = gfp_stack;
end

%%
index = 10;
trange =103:143;
ct = min(trange)+5;
gfp_stack = gpf_cell{index};
cimage = gfp_stack(:,:,ct);
cimage = cimage./2000;
[cropeedimage,rectout] = imcrop(cimage);
close all
%%
removesmallarea = 30;
removebigarea = 500;
lowintensitylimit = 000;
ratio = 0.5;

cubeout = images.spatialref.Cuboid([rectout(1),rectout(1)+rectout(3)],...
                                   [rectout(2),rectout(2)+rectout(4)],...
                                   [trange(1),trange(end)]);
stackofinterest = imcrop3(gfp_stack,cubeout);
intensitylist = [];
montagecheckimage = [];
for ct = trange-trange(1)+1
    %Segmentation
    cimage = stackofinterest(:,:,ct);
    cmax = max(cimage,[],'all');
    cmax(cmax < lowintensitylimit) = lowintensitylimit;
    cmask = cimage > ratio*cmax;
    csmallmask = cmask-bwareaopen(cmask,removesmallarea) == 1;
    cimage(csmallmask) = 0;
    cmax = max(cimage,[],'all');
    cmax(cmax < lowintensitylimit) = lowintensitylimit;
    cmask = cimage > ratio*cmax;
    cmask = bwareaopen(cmask,removesmallarea,4);
    cmask(bwareaopen(cmask,removebigarea)) = 0;

    Label = labelmatrix(bwconncomp(cmask));
    Nucleiecc = regionprops(cmask,"Circularity");
    thinnuclei = find([Nucleiecc.Circularity] < 0.3);
    for i = thinnuclei
        Label(Label==i) = 0;
    end
    cmask=cmask&Label;
%    [~,maxspot] = max(cimage,[],'all');
%    [maxint_x,maxint_y] = ind2sub(size(cimage),maxspot);
%    boundary = bwtraceboundary(cmask,[maxint_x,maxint_y],'N');
%     hold on;
%     imshow(cimage./2000)
%     plot(boundary(:,2),boundary(:,1),'g','LineWidth',3);
%     hold off

    %montage generation
    cimage = stackofinterest(:,:,ct);
    cmontageimage = cat(3,cimage,cimage,cimage);
    ccolormask = cat(3,cimage<0,cmask,cimage<0);
    cmontageimage(ccolormask) = 4000;
    montagecheckimage{ct} = cmontageimage./4000;

    %Intensity
    intensitylist = [intensitylist mean(cimage(cmask))];
end
intensitylist(isnan(intensitylist))=0;

close all
figure,montage(montagecheckimage)
figure,plot(intensitylist)
xlabel('Time (Hours)')
ylabel('Raw Magnitude')

%%

intensitycellsfilename = [segdir '\seg.mat'];
intindex = 1;
intensitycells = [];
if exist(intensitycellsfilename,'file')
    load(intensitycellsfilename) 
    intindex = length(intensitycells)+1;
end
if exist('intensitylist','var')
    intensitycells{intindex} = intensitylist;
    save(intensitycellsfilename,'intensitycells')
    save([segdir '\segindex' num2str(intindex) 'Slide' num2str(index) 'Loc' num2str(rectout(1)) '_' num2str(rectout(2)) '.mat'])
end

%%
%raw plot
basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF';
savedir = [basedir '\Analysis\gfp'];
segdir = [basedir '\Analysis\seg'];
celltypefile = [segdir '\celltypeinfo.mat'];
if exist(celltypefile,'file')
    load(celltypefile)
end
intensitycellsfilename = [segdir '\seg.mat'];
if exist(intensitycellsfilename,'file')
    load(intensitycellsfilename) 
end
arealistfile = [segdir '\arealistinfo.mat'];
if exist(arealistfile,'file')
    load(arealistfile) 
end

clf
figure(1)
hold on
for i = 1:length(intensitycells)
    cint = intensitycells{i};
    cevent = find(cint==0, 1, 'last' )-1;
    x = ((0:length(cint)-1)-cevent)*10;
    y = double(cint);
    if strcmp(celltype{i},'Primary')
        color = 'red';
    else 
        color = 'blue';
    end
    plot(x/60,y,color)
end
hold off
xlabel('Time (Hours)')
ylabel('Raw Magnitude')
%normalized

figure(2)
hold on
endt = [];
maxintensity = [];
recorded = false;
for i = 1:length(intensitycells)
    cint = intensitycells{i};
    cint(cint==0) = [];
    x = 0:10:length(cint)*10-10;
    y = double(cint./max(cint));

    belowcutoff = y<0.7;
    [maxin,maxindex] = max(y);
    maxintensity(i) = max(cint);
    belowcutoff(1:maxindex) = 0;
    index = find(belowcutoff,1);
    endt = [endt x(index)];
    if strcmp(celltype{i},'Primary')
        color = 'red';
%         plot(x/60,y,color)
%         xlim([0 18])
    else 
        color = 'blue';
%         plot(x/60,y,color)
%         xlim([0 18])  
    end
    plot(x/60,y,color)
end
hold off
xlabel('Time (Hours)')
ylabel('Normalized Magnitude')
xlim([0 22]);
title('Tracking Rupture BAF Intensity Over Time')
legend('Micronucleus','Primary Nucleus');
exportgraphics(gca,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Analysis\Graphs\Tracking Rupture BAF Intensity Over Time.pdf'],'Resolution',600,'ContentType','vector')

figure(3)
hold on
endt = [];
recorded = false;
arealistjet = round(arealist./max(arealist)*256);
maxinlistjet = round(maxintensity./max(maxintensity)*256);
for i = 1:length(intensitycells)
    cint = intensitycells{i};
    cint(cint==0) = [];
    x = 0:10:length(cint)*10-10;
    y = double(cint./max(cint));

    belowcutoff = y<0.7;
    [maxin,maxindex] = max(y);
    belowcutoff(1:maxindex) = 0;
    index = find(belowcutoff,1);
    endt = [endt x(index)];
    jet = jet;
    color = jet(maxinlistjet(i),:);
    plot(x/60,y,'Color',color)
end
hold off
xlabel('Time (Hours)')
ylabel('Normalized Magnitude')
colorbar('TickLabels','')
colormap("jet")
exportgraphics(gca,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Analysis\Graphs\Tracking Rupture BAF Intensity Heatmap Intensity.pdf'],'Resolution',600,'ContentType','vector')

figure(31)
hold on
endt = [];
recorded = false;
arealistjet = round(arealist./max(arealist)*256);
maxinlistjet = round(maxintensity./max(maxintensity)*256);
for i = 1:length(intensitycells)
    cint = intensitycells{i};
    cint(cint==0) = [];
    x = 0:10:length(cint)*10-10;
    y = double(cint./max(cint));

    belowcutoff = y<0.7;
    [maxin,maxindex] = max(y);
    belowcutoff(1:maxindex) = 0;
    index = find(belowcutoff,1);
    endt = [endt x(index)];
    jet = jet;
    color = jet(arealistjet(i),:);
    plot(x/60,y,'Color',color)
end
hold off
xlabel('Time (Hours)')
ylabel('Normalized Magnitude')
colorbar('TickLabels','')
colormap("jet")
exportgraphics(gca,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Analysis\Graphs\Tracking Rupture BAF Intensity Heatmap Area.pdf'],'Resolution',600,'ContentType','vector')


[~,sortedindex] = sort(endt);
a = [num2cell(sortedindex)',num2cell(endt(sortedindex))',celltype(sortedindex)',num2cell(arealist(sortedindex))',num2cell(maxintensity(sortedindex))'];
a = [{'Index','Time','Classification','Area','Max Intensity'};a];

figure(4)
hold on
endt = [];
maxintensity = [];
recorded = false;
for i = 1:length(intensitycells)
    cint = intensitycells{i};
    cint(cint==0) = [];
    x = 0:10:length(cint)*10-10;
    y = double(cint./max(cint));

    belowcutoff = y<0.7;
    [maxin,maxindex] = max(y);
    maxintensity(i) = max(cint);
    belowcutoff(1:maxindex) = 0;
    index = find(belowcutoff,1);
    endt = [endt x(index)];
    if strcmp(celltype{i},'Primary')
        color = 'red';
%         plot(x/60,y,color)
%         xlim([0 18])
    else 
        color = 'blue';
        plot(x/60,y,color)
        xlim([0 18])  
    end
%    plot(x/60,y,color)
end
hold off
xlabel('Time (Hours)')
ylabel('Normalized Magnitude')
xlim([0 22]);
title('Tracking Rupture BAF Intensity Over Time')
legend('Micronucleus');
exportgraphics(gca,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Analysis\Graphs\Tracking Rupture BAF Intensity Over Time MNonly.pdf'],'Resolution',600,'ContentType','vector')

figure(5)
hold on
endt = [];
maxintensity = [];
recorded = false;
for i = 1:length(intensitycells)
    cint = intensitycells{i};
    cint(cint==0) = [];
    x = 0:10:length(cint)*10-10;
    y = double(cint./max(cint));

    belowcutoff = y<0.7;
    [maxin,maxindex] = max(y);
    maxintensity(i) = max(cint);
    belowcutoff(1:maxindex) = 0;
    index = find(belowcutoff,1);
    endt = [endt x(index)];
    if strcmp(celltype{i},'Primary')
        color = 'red';
        plot(x/60,y,color)
        xlim([0 18])
    else 
        color = 'blue';
%         plot(x/60,y,color)
%         xlim([0 18])  
    end
%    plot(x/60,y,color)
end
hold off
xlabel('Time (Hours)')
ylabel('Normalized Magnitude')
xlim([0 22]);
title('Tracking Rupture BAF Intensity Over Time')
legend('Primary Nucleus');
exportgraphics(gca,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Analysis\Graphs\Tracking Rupture BAF Intensity Over Time PNonly.pdf'],'Resolution',600,'ContentType','vector')

%% Graphing
figure(4),
hold on;
celltypemicroistrue = [];
for i = 1:length(celltype)
    celltypemicroistrue(i) = strcmp(celltype{i},'Micro');
end
eventtime = [a{2:end,2}]./60;
arealist = [a{2:end,4}];
x = eventtime;
y = arealist;
% Get coefficients of a line fit through the data.
[core,prob] = corrcoef(eventtime,arealist)
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
for choice = 0:1
    ceventtime = eventtime(celltypemicroistrue == choice);
    carealist = arealist(celltypemicroistrue == choice);
    scatter(ceventtime,carealist,"filled")
    %corrcoef(ceventtime,carealist)
    xlabel('Time (Hours)')
    ylabel('Event Area (Pixels)')
end
legend('Line of Best Fit','Micronucleus','Primary Nucleus')
legend boxoff   
ylim([0 max(carealist)*1.3])
hold off;
exportgraphics(gca,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Analysis\Graphs\LineofBFArea.pdf'],'Resolution',600,'ContentType','vector')

figure(5),
hold on;
celltypemicroistrue = [];
for i = 1:length(celltype)
    celltypemicroistrue(i) = strcmp(celltype{i},'Micro');
end
eventtime = [a{2:end,2}]./60;
maxinlist = [a{2:end,5}];
x = eventtime;
y = maxinlist;
% Get coefficients of a line fit through the data.
[core,prob] = corrcoef(eventtime,maxinlist)
coefficients = polyfit(x, y, 1);
xFit = linspace(min(x), max(x), 1000);
yFit = polyval(coefficients , xFit);
plot(xFit, yFit, 'r-', 'LineWidth', 2); % Plot fitted line.
for choice = 0:1
    ceventtime = eventtime(celltypemicroistrue == choice);
    cmaxinlist = maxinlist(celltypemicroistrue == choice);
    scatter(ceventtime,cmaxinlist,"filled")
    %corrcoef(ceventtime,carealist)
    xlabel('Time (Hours)')
    ylabel('Max Intensity')
end
legend('Line of Best Fit','Micronucleus','Primary Nucleus')
legend boxoff   
ylim([0 max(cmaxinlist)*1.3])
hold off;
exportgraphics(gca,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Analysis\Graphs\LineofBFInten.pdf'],'Resolution',600,'ContentType','vector')

%% Output Data Table
basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF';
savedir = [basedir '\Analysis\gfp'];
segdir = [basedir '\Analysis\seg'];
celltypefile = [segdir '\celltypeinfo.mat'];
if exist(celltypefile,'file')
    load(celltypefile)
end
intensitycellsfilename = [segdir '\seg.mat'];
if exist(intensitycellsfilename,'file')
    load(intensitycellsfilename) 
end
arealistfile = [segdir '\arealistinfo.mat'];
if exist(arealistfile,'file')
    load(arealistfile) 
end
filenamelist = dir(segdir);
indexlist = [];
Slidelist = [];
loc_xlist = [];
loc_ylist = [];
time_start = [];
for i = 1:length(filenamelist)
    cfilename = filenamelist(i).name;
    if contains(cfilename,'segindex')
        [temp,rest] = strtok(cfilename,'S');
        cindex = str2num(temp(9:end));
        indexlist = [indexlist cindex];
        [temp,rest] = strtok(rest,'L');
        cSlide = str2num(temp(6:end));
        Slidelist = [Slidelist cSlide];
        [temp,rest] = strtok(rest,'_');
        cloc_x = str2num(temp(4:end));
        loc_xlist = [loc_xlist cloc_x];
        [temp,rest] = strtok(rest,'m');
        cloc_y = str2num(temp(2:end-1));
        loc_ylist = [loc_ylist cloc_y];

        load([segdir '\' cfilename],'cubeout')
        time_start = [time_start cubeout.ZLimits(1)];
    end
end
endt = [];
maxintensity = [];
meanintensity = [];
for i = 1:length(intensitycells)
    cint = intensitycells{i};
    cint(cint==0) = [];
    x = 0:10:length(cint)*10-10;
    y = double(cint./max(cint));

    belowcutoff = y<0.7;
    [maxin,maxindex] = max(y);
    maxintensity(i) = max(cint);
    meanintensity(i) = mean(cint);
    belowcutoff(1:maxindex) = 0;
    index = find(belowcutoff,1);
    endt = [endt x(index)];
end
time_start = time_start * 10;
time_end = time_start + endt;
sortedindex = 1:length(endt);
a = [num2cell(sortedindex)',num2cell(Slidelist)',num2cell(endt(sortedindex))',num2cell(time_start)',num2cell(time_end)'...
    celltype(sortedindex)',num2cell(arealist(sortedindex))',...
    num2cell(maxintensity(sortedindex))',num2cell(meanintensity(sortedindex))',num2cell(loc_xlist)',num2cell(loc_ylist)'];
T = cell2table(a,"VariableNames",...
    {'Index','SlideNumber','Time(mins)','TimeStart','TimeEnd'...
    'Classification','Area',...
    'Max BAF Intensity','Mean BAF Intensity','X_Coordinates','Y_Coordinates'});
writetable(T,'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF\Analysis\seg\Kuramochi20180524RuptureTrackingMetaData.xlsx');
%% N and p between MN PN

PN_time = [];
MN_time = [];
endt = [];
maxintensity = [];
for i = 1:length(intensitycells)
    cint = intensitycells{i};
    cint(cint==0) = [];
    x = 0:10:length(cint)*10-10;
    y = double(cint./max(cint));
    belowcutoff = y<0.7;
    [maxin,maxindex] = max(y);
    maxintensity(i) = max(cint);
    belowcutoff(1:maxindex) = 0;
    index = find(belowcutoff,1);
    endt = [endt x(index)];
    if strcmp(celltype{i},'Primary')
        PN_time = [PN_time x(index)];
    elseif strcmp(celltype{i},'Micro')
        MN_time = [MN_time x(index)];
    end
end

mean(PN_time)
mean(MN_time)
[~,p] = ttest2(PN_time,MN_time,'Tail','left');

%% Individual review
basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF';
savedir = [basedir '\Analysis\gfp'];
segdir = [basedir '\Analysis\seg'];

segnamelist = {dir(segdir).name};
segnamelist = segnamelist(contains(segnamelist,'index'));
for namei = 1:length(segnamelist)
    csegname = segnamelist{namei};
    [num,~] = strtok(csegname(9:end),'Slide');
    indexsegmaskwhatever(namei) = str2num(num);
end
[~,sindex] = sort(indexsegmaskwhatever);
segnamelist = segnamelist(sindex);
lastindex = 0;
celltype = {};
for namei = 1:length(segnamelist)
    csegname = segnamelist{namei};
    [num,~] = strtok(csegname(9:end),'Slide');
    num = str2num(num);
    
    load([segdir '\' csegname])
    celltypefile = [segdir '\celltypeinfo.mat'];
    if exist(celltypefile,'file')
        load(celltypefile)
    end
    if lastindex ~= index
        gfp_stack = gpf_cell{index};
        lastindex = index;
    end

    ct = min(trange)+5;
    cimage = gfp_stack(:,:,ct);
    cimage = cimage./2000;
    hold on
    imshow(cimage)
    rectangle('Position',rectout,'LineWidth',2,'Edgecolor','green')
    hold off
    ct 
    index
    answer = questdlg('Primary Nucleus?', ...
	'Selection', ...
	'Primary','Micro','cancel');
    celltype{num} = answer;
    %save(celltypefile, 'celltype')
    close all
end

%% Area Collection
basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF';
savedir = [basedir '\Analysis\gfp'];
segdir = [basedir '\Analysis\seg'];

segnamelist = {dir(segdir).name};
segnamelist = segnamelist(contains(segnamelist,'index'));
arealist = [];
intensitycells_new = [];
indexsegmaskwhatever = [];
for namei = 1:length(segnamelist)
    csegname = segnamelist{namei};
    [num,~] = strtok(csegname(9:end),'Slide');
    indexsegmaskwhatever(namei) = str2num(num);
end
[~,sindex] = sort(indexsegmaskwhatever);
segnamelist = segnamelist(sindex);

for namei = 1:length(segnamelist)
    csegname = segnamelist{namei};
    load([segdir '\' csegname])
    carealist = [];
    csegname = segnamelist{namei};
    [num,~] = strtok(csegname(9:end),'Slide');
    num = str2num(num)
    ct = min(trange)+5;
    gfp_stack = gpf_cell{index};
    if namei <= 2
        %some early images need to be transpose due to an error
         gfp_stack=pagetranspose( gfp_stack);
    end
    
    cubeout = images.spatialref.Cuboid([rectout(1),rectout(1)+rectout(3)],...
                                       [rectout(2),rectout(2)+rectout(4)],...
                                       [trange(1),trange(end)]);
    stackofinterest = imcrop3(gfp_stack,cubeout);
    intensitylist = [];
    montagecheckimage = [];
    if ~exist('ratio','var')
        ratio = 0.4;
    end
    for ct = trange-trange(1)+1
        %Segmentation
        cimage = stackofinterest(:,:,ct);
        cmax = max(cimage,[],'all');
        cmax(cmax < lowintensitylimit) = lowintensitylimit;
        cmask = cimage > ratio*cmax;
        csmallmask = cmask-bwareaopen(cmask,removesmallarea) == 1;
        cimage(csmallmask) = 0;
        cmax = max(cimage,[],'all');
        cmax(cmax < lowintensitylimit) = lowintensitylimit;
        cmask = cimage > ratio*cmax;
        cmask = bwareaopen(cmask,removesmallarea,4);
        cmask(bwareaopen(cmask,removebigarea)) = 0;
    
        Label = labelmatrix(bwconncomp(cmask));
        Nucleiecc = regionprops(cmask,"Circularity");
        thinnuclei = find([Nucleiecc.Circularity] < 0.2);
        for i = thinnuclei
            Label(Label==i) = 0;
        end
        cmask=cmask&Label;
    %    [~,maxspot] = max(cimage,[],'all');
    %    [maxint_x,maxint_y] = ind2sub(size(cimage),maxspot);
    %    boundary = bwtraceboundary(cmask,[maxint_x,maxint_y],'N');
    %     hold on;
    %     imshow(cimage./2000)
    %     plot(boundary(:,2),boundary(:,1),'g','LineWidth',3);
    %     hold off
        carealist = [carealist sum(cmask,'all')];
        %montage generation
        cimage = stackofinterest(:,:,ct);
        cimage = cimage*3;
        cimage = cimage-mean(cimage)*0;
        cimage(cimage<0) = 0;
        cmontageimage = cat(3,cimage,cimage,cimage);
        ccolormask = cat(3,cimage<0,cmask,cimage<0);
        cmontageimage(ccolormask) = 4000;
        montagecheckimage{ct} = cmontageimage./4000;
    
        %Intensity
        intensitylist = [intensitylist mean(cimage(cmask))];
    end
    intensitylist(isnan(intensitylist))=0;
    intensitycells_new{num} = intensitylist;
    figure,montage(montagecheckimage)
    figure,plot(intensitylist)
    xlabel('Time (Hours)')
    ylabel('Raw Magnitude')

    carealist(carealist == 0) = [];
    carealist = carealist(1:10);
    careamean = round(mean(carealist));
    arealist(num) = careamean;
end
intensitycells = intensitycells_new;
close all
arealistfile = [segdir '\arealistinfo.mat'];
%save(arealistfile, 'arealist')


%% Sample generation
basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\COY_Old_LiveCell\Live_Cell_Vids_from_Alex\18_0524 Kuramochi GFP-BAF';
savedir = [basedir '\Analysis\gfp'];
segdir = [basedir '\Analysis\seg'];

segnamelist = {dir(segdir).name};
segnamelist = segnamelist(contains(segnamelist,'index'));
arealist = [];
intensitycells_new = [];
indexsegmaskwhatever = [];
for namei = 1:length(segnamelist)
    csegname = segnamelist{namei};
    [num,~] = strtok(csegname(9:end),'Slide');
    indexsegmaskwhatever(namei) = str2num(num);
end
[~,sindex] = sort(indexsegmaskwhatever);
segnamelist = segnamelist(sindex);

for namei = 3
    csegname = segnamelist{namei};
    load([segdir '\' csegname])
    carealist = [];
    csegname = segnamelist{namei};
    [num,~] = strtok(csegname(9:end),'Slide');
    num = str2num(num)
    ct = min(trange)+5;
    gfp_stack = gpf_cell{index};
    if namei <= 2
        %some early images need to be transpose due to an error
         gfp_stack=pagetranspose( gfp_stack);
    end
    
    cubeout = images.spatialref.Cuboid([rectout(1),rectout(1)+rectout(3)],...
                                       [rectout(2),rectout(2)+rectout(4)],...
                                       [trange(1),trange(end)]);
    stackofinterest = imcrop3(gfp_stack,cubeout);
    intensitylist = [];
    montagecheckimage = [];
    if ~exist('ratio','var')
        ratio = 0.4;
    end
    for ct = trange-trange(1)+1
        %Segmentation
        cimage = stackofinterest(:,:,ct);
        cmax = max(cimage,[],'all');
        cmax(cmax < lowintensitylimit) = lowintensitylimit;
        cmask = cimage > ratio*cmax;
        csmallmask = cmask-bwareaopen(cmask,removesmallarea) == 1;
        cimage(csmallmask) = 0;
        cmax = max(cimage,[],'all');
        cmax(cmax < lowintensitylimit) = lowintensitylimit;
        cmask = cimage > ratio*cmax;
        cmask = bwareaopen(cmask,removesmallarea,4);
        cmask(bwareaopen(cmask,removebigarea)) = 0;
    
        Label = labelmatrix(bwconncomp(cmask));
        Nucleiecc = regionprops(cmask,"Circularity");
        thinnuclei = find([Nucleiecc.Circularity] < 0.2);
        for i = thinnuclei
            Label(Label==i) = 0;
        end
        cmask=cmask&Label;
    %    [~,maxspot] = max(cimage,[],'all');
    %    [maxint_x,maxint_y] = ind2sub(size(cimage),maxspot);
    %    boundary = bwtraceboundary(cmask,[maxint_x,maxint_y],'N');
    %     hold on;
    %     imshow(cimage./2000)
    %     plot(boundary(:,2),boundary(:,1),'g','LineWidth',3);
    %     hold off
        carealist = [carealist sum(cmask,'all')];
        %montage generation
        cimage = stackofinterest(:,:,ct);
        cimage = cimage*3;
        cimage = cimage-mean(cimage,'all')*0.5;
        cimage(cimage<0) = 0;
        cmontageimage = cat(3,cimage,cimage,cimage);
        ccolormask = cat(3,cimage<0,cmask,cimage<0);
        cmontageimage(ccolormask) = 4000;
        montagecheckimage{ct} = cmontageimage./4000;
    
        %Intensity
        intensitylist = [intensitylist mean(cimage(cmask))];
    end
    intensitylist(isnan(intensitylist))=0;
    intensitycells_new{num} = intensitylist;
    figure,montage(montagecheckimage)
    figure,plot(intensitylist)
    xlabel('Time (Hours)')
    ylabel('Raw Magnitude')

    carealist(carealist == 0) = [];
    carealist = carealist(1:10);
    careamean = round(mean(carealist));
    arealist(num) = careamean;
end






