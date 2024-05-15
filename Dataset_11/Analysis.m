%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primary Nucleus Segmentation (Morphology Optimized)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\';
csvoutputdir = [basedir 'CSV_output\'];
csvoutputdatasetdir = [csvoutputdir 'GBM_84_Region1\'];
filenamestart = 'GBM_Registered_FC3';
csvfilelist = [csvoutputdatasetdir filenamestart '_Statistics\'];

Results_morphology = [];

markers = {'Hoechst1','BANF1','STING','TREX1'...
           'Hoechst2','NFkB','pH2AX','IRF3'...
           'Hoechst3','SOX2','cGAS','SOX9'...
           'Hoechst4','ELAVL4','DCX','OLIG2'...
           'Hoechst5','EGFR','CEBPB','SOX10'...
           'Hoechst6','TMEM119','LAMINAC','HOPX'...
           'Hoechst7','Ki67','LaminB1','LaminB2'...
           'Hoechst8','HIF1a','NDRG1','IFITM1'...
           'Hoechst9','CD31','MX1','INSM1'...
           'Hoechst10','BANF1_2','cGAS_2','pH2AX_2'...
            }; 
Lineages = {'SOX9','ELAVL4','HIF1a','DCX','OLIG2','SOX10'};
Interferons = {'NFkB','MX1','IRF3','IFITM1'...
            };
% Channel Mean
for i = 1:40
    cmeanfile = [csvfilelist filenamestart '_Intensity_Mean_Ch=' num2str(i) '_Img=1.csv'];
    ctable = readtable(cmeanfile);
    markersname = markers{i};
    markersname(markersname == 'κ') = 'k';
    for ri = 1:length(ctable.IntensityMean)
        Results_morphology(ri).(markersname) = ctable.IntensityMean(ri);
    end
end

% Channel std
for i = 1:40
    cmeanfile = [csvfilelist filenamestart '_Intensity_StdDev_Ch=' num2str(i) '_Img=1.csv'];
    ctable = readtable(cmeanfile);
    markersname = markers{i};
    markersname(markersname == 'κ') = 'k';
    markersname = [markersname '_std'];
    for ri = 1:length(ctable.IntensityStdDev)
        Results_morphology(ri).(markersname) = ctable.IntensityStdDev(ri);
    end
end

filelistall = dir(csvfilelist);
filenamelist = {filelistall.name};
for i = 1:length(filenamelist)
    if contains(filenamelist{i},'Sphericity.csv') || ...
            contains(filenamelist{i},'Ellipticity_(oblate).csv') || ...
            contains(filenamelist{i},'Ellipticity_(prolate).csv') || ...
            contains(filenamelist{i},'Volume.csv') || ...
            contains(filenamelist{i},'Area.csv') || ...
            contains(filenamelist{i},'StdDev_Ch=1_Img=1.csv') 
        csvoutcell = readcell([csvfilelist filenamelist{i}]);
        [mrow,~] = size(csvoutcell);
        cfieldname = csvoutcell{1,1};
        cfieldname = regexprep(cfieldname, {'[%() ]+', '_+$'}, {'', ''});
        if strcmp(cfieldname,'Intensity')
            cfieldname = 'StdDev';
        end
        for rowi = 4:mrow
            if ischar(csvoutcell{rowi,1})
                Results_morphology(rowi-3).(cfieldname) = str2num(strtok(csvoutcell{rowi,1},','));
            else
                Results_morphology(rowi-3).(cfieldname) = csvoutcell{rowi,1};
            end
            %Results(rowi-3+appending_index).BAF = bafdiri == 2;
        end
    end
end

% Morphology Score
logVolume = log10([Results_morphology.Volume]);
nVolume = (logVolume-prctile(logVolume,5))./(prctile(logVolume,95)-prctile(logVolume,5));
%nVolume  = [Results_morphology.Volume];
IrregularIndex = nVolume ./[Results_morphology.Sphericity];
for ri = 1:length([Results_morphology.Volume])
    Results_morphology(ri).IrregularIndex = IrregularIndex(ri);
end


% % Remove error cells
% deletelist = [];
% for ri = 1:length([Results_morphology.Volume])
%     if Results_morphology(ri).Volume < 1
%         deletelist = [deletelist ri];
%     end
% end
% Results_morphology(deletelist) = [];
% BAF Pos vs Neg
bafmaxfile = [csvfilelist filenamestart '_Intensity_Max_Ch=2_Img=1'];
csvoutcell = readtable(bafmaxfile);
bafposmask = csvoutcell.IntensityMax > 3160;
for ri = 1:length([Results_morphology.Volume])
    Results_morphology(ri).Rupture = bafposmask(ri);
end
% SOX2 Mask
sox2mask = [Results_morphology.SOX2]' > 300;
for ri = 1:length([Results_morphology.Volume])
    Results_morphology(ri).SOX2mask = sox2mask(ri);
end
% Location
for axis = {'X','Y','Z'}
    locfile = [csvfilelist filenamestart '_Position_' axis{1} '.csv'];
    csvoutcell = readtable(locfile);
    cpos = table2array([csvoutcell(:,1)]);
    for ri = 1:length([Results_morphology.Volume])
        Results_morphology(ri).(axis{1}) = cpos(ri);
    end
end

%Save Results
saveresult = struct2table(Results_morphology);
writetable(saveresult,[basedir '3D_Morphology.csv'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Micronucleus Segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\';
csvoutputdir = [basedir 'CSV_output\'];
csvoutputdatasetdir = [csvoutputdir 'GBM_84_Region1\'];
filenamestart = 'GBM_MNMN_Registered';
csvfilelist = [csvoutputdatasetdir filenamestart '_Statistics\'];

Results_micro = [];

% Channel Mean
for i = 1:40
    cmeanfile = [csvfilelist filenamestart '_Intensity_Mean_Ch=' num2str(i) '_Img=1.csv'];
    ctable = readtable(cmeanfile);
    markersname = markers{i};
    markersname(markersname == 'κ') = 'k';
    for ri = 1:length(ctable.IntensityMean)
        Results_micro(ri).(markersname) = ctable.IntensityMean(ri);
    end
end

appending_index = 0;

filelistall = dir(csvfilelist);
filenamelist = {filelistall.name};
for i = 1:length(filenamelist)
    if contains(filenamelist{i},'Sphericity.csv') || ...
            contains(filenamelist{i},'Ellipticity_(oblate).csv') || ...
            contains(filenamelist{i},'Ellipticity_(prolate).csv') || ...
            contains(filenamelist{i},'Volume.csv') || ...
            contains(filenamelist{i},'Area.csv') 
        csvoutcell = readcell([csvfilelist filenamelist{i}]);
        [mrow,~] = size(csvoutcell);
        cfieldname = csvoutcell{1,1};
        cfieldname = regexprep(cfieldname, {'[%() ]+', '_+$'}, {'', ''});
        for rowi = 4:mrow
            if ischar(csvoutcell{rowi,1})
                Results_micro(rowi-3).(cfieldname) = str2num(strtok(csvoutcell{rowi,1},','));
            else
                Results_micro(rowi-3).(cfieldname) = csvoutcell{rowi,1};
            end
            %Results(rowi-3+appending_index).BAF = bafdiri == 2;
        end
    end
end

for axis = {'X','Y','Z'}
    locfile = [csvfilelist filenamestart '_Position_' axis{1} '.csv'];
    csvoutcell = readtable(locfile);
    cpos = table2array([csvoutcell(:,1)]);
    for ri = 1:length([Results_micro.Volume])
        Results_micro(ri).(axis{1}) = cpos(ri);
    end
end

bafposmask = [Results_micro.BANF1]' > 250;
for ri = 1:length([Results_micro.Volume])
    Results_micro(ri).Rupture = bafposmask(ri);
end
% % SOX2 Mask
% sox2mask = [Results_micro.SOX2]' > 30;
% for ri = 1:length([Results_micro.Volume])
%     Results_micro(ri).SOX2mask = sox2mask(ri);
% end

% Save Results
saveresult = struct2table(Results_micro);
writetable(saveresult,[basedir '3D_MN.csv'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Primary Nucleus Segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\';
csvoutputdir = [basedir 'CSV_output\'];
csvoutputdatasetdir = [csvoutputdir 'GBM_84_Region1\'];
filenamestart = 'GBM_Registered_All_Statistics';
csvfilelist = [csvoutputdatasetdir filenamestart '_Statistics\'];

Results_all = [];

% Channel Mean
for i = 1:40
    cmeanfile = [csvfilelist filenamestart '_Intensity_Mean_Ch=' num2str(i) '_Img=1.csv'];
    ctable = readtable(cmeanfile);
    markersname = markers{i};
    markersname(markersname == 'κ') = 'k';
    for ri = 1:length(ctable.IntensityMean)
        Results_all(ri).(markersname) = ctable.IntensityMean(ri);
    end
end

filelistall = dir(csvfilelist);
filenamelist = {filelistall.name};
for i = 1:length(filenamelist)
    if contains(filenamelist{i},'Sphericity.csv') || ...
            contains(filenamelist{i},'Ellipticity_(oblate).csv') || ...
            contains(filenamelist{i},'Ellipticity_(prolate).csv') || ...
            contains(filenamelist{i},'Volume.csv') || ...
            contains(filenamelist{i},'Area.csv') 
        csvoutcell = readcell([csvfilelist filenamelist{i}]);
        [mrow,~] = size(csvoutcell);
        cfieldname = csvoutcell{1,1};
        cfieldname = regexprep(cfieldname, {'[%() ]+', '_+$'}, {'', ''});
        for rowi = 4:mrow
            if ischar(csvoutcell{rowi,1})
                Results_all(rowi-3).(cfieldname) = str2num(strtok(csvoutcell{rowi,1},','));
            else
                Results_all(rowi-3).(cfieldname) = csvoutcell{rowi,1};
            end
            %Results(rowi-3+appending_index).BAF = bafdiri == 2;
        end
    end
end

% Morphology Score
logVolume = log10([Results_all.Volume]);
nVolume = (logVolume-prctile(logVolume,5))./prctile(logVolume,95);
%nVolume = ([Results_all.Volume]-prctile([Results_all.Volume],0))./prctile([Results_all.Volume],95);
%nVolume = log([Results_all.Volume]+1);
%nVolume  = [Results_all.Volume];
IrregularIndex = nVolume ./[Results_all.Sphericity] ;
for ri = 1:length([Results_all.Volume])
    Results_all(ri).IrregularIndex = IrregularIndex(ri);
end

% Location
for axis = {'X','Y','Z'}
    locfile = [csvfilelist filenamestart '_Position_' axis{1} '.csv'];
    csvoutcell = readtable(locfile);
    cpos = table2array([csvoutcell(:,1)]);
    for ri = 1:length(Results_all)
        Results_all(ri).(axis{1}) = cpos(ri);
    end
end

emptymask = find([Results_all.Hoechst6] == 0);
Results_all(emptymask) = [];

% SOX2 Mask
sox2mask = [Results_all.SOX2]' > 200;
for ri = 1:length(Results_all)
    Results_all(ri).SOX2mask = sox2mask(ri);
end

%Save Results
saveresult = struct2table(Results_all);
writetable(saveresult,[basedir '3D_all.csv'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAF Foci Segmentation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

basedir = 'Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\';
csvoutputdir = [basedir 'CSV_output\'];
csvoutputdatasetdir = [csvoutputdir 'GBM_84_Region1\'];
filenamestart = 'GBM_Registered_AllCells_BAFFoci';
csvfilelist = [csvoutputdatasetdir filenamestart '_Statistics\'];

Results_Foci = [];

% Channel Mean
for i = 1:40
    cmeanfile = [csvfilelist filenamestart '_Intensity_Mean_Ch=' num2str(i) '_Img=1.csv'];
    ctable = readtable(cmeanfile);
    markersname = markers{i};
    markersname(markersname == 'κ') = 'k';
    for ri = 1:length(ctable.IntensityMean)
        Results_Foci(ri).(markersname) = ctable.IntensityMean(ri);
    end
end

% Location
for axis = {'X','Y','Z'}
    locfile = [csvfilelist filenamestart '_Position_' axis{1} '.csv'];
    csvoutcell = readtable(locfile);
    cpos = table2array([csvoutcell(:,1)]);
    for ri = 1:length(Results_Foci)
        Results_Foci(ri).(axis{1}) = cpos(ri);
    end
end

appending_index = 0;

filelistall = dir(csvfilelist);
filenamelist = {filelistall.name};
for i = 1:length(filenamelist)
    if contains(filenamelist{i},'Sphericity.csv') || ...
            contains(filenamelist{i},'Ellipticity_(oblate).csv') || ...
            contains(filenamelist{i},'Ellipticity_(prolate).csv') || ...
            contains(filenamelist{i},'Volume.csv') || ...
            contains(filenamelist{i},'Area.csv') 
        csvoutcell = readcell([csvfilelist filenamelist{i}]);
        [mrow,~] = size(csvoutcell);
        cfieldname = csvoutcell{1,1};
        cfieldname = regexprep(cfieldname, {'[%() ]+', '_+$'}, {'', ''});
        for rowi = 4:mrow
            if ischar(csvoutcell{rowi,1})
                Results_Foci(rowi-3).(cfieldname) = str2num(strtok(csvoutcell{rowi,1},','));
            else
                Results_Foci(rowi-3).(cfieldname) = csvoutcell{rowi,1};
            end
            %Results(rowi-3+appending_index).BAF = bafdiri == 2;
        end
    end
end

% Save Results
saveresult = struct2table(Results_Foci);
writetable(saveresult,[basedir '3D_Foci.csv'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BAF Foci, MN mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Map Foci to Primary nucleus
xlist = [Results_all.X];
ylist = [Results_all.Y];
zlist = [Results_all.Z];
coordlist = [xlist;ylist;zlist]';
xlist_Foci = [Results_Foci.X];
ylist_Foci = [Results_Foci.Y];
zlist_Foci = [Results_Foci.Z];
coordlist_Foci = [xlist_Foci;ylist_Foci;zlist_Foci]';
[k_Foci,dist] = dsearchn(coordlist,coordlist_Foci);
for i = 1:length(Results_all)
    fociindexlist = find(k_Foci==i);
    Results_all(i).FociN = length(fociindexlist);
end
% % Map Foci to Primary nucleus
% xlist = [Results_morphology.X];
% ylist = [Results_morphology.Y];
% zlist = [Results_morphology.Z];
% coordlist = [xlist;ylist;zlist]';
% xlist_Foci = [Results_Foci.X];
% ylist_Foci = [Results_Foci.Y];
% zlist_Foci = [Results_Foci.Z];
% coordlist_Foci = [xlist_Foci;ylist_Foci;zlist_Foci]';
% [k_Foci,dist] = dsearchn(coordlist,coordlist_Foci);
% for i = 1:length(Results_morphology)
%     toofarcount = sum(dist(k_Foci==i)>10);
%     fociindexlist = find(k_Foci==i);
%     Results_morphology(i).FociN = length(fociindexlist)-toofarcount;
% end

% Map MN to Primary nucleus
xlist = [Results_all.X];
ylist = [Results_all.Y];
zlist = [Results_all.Z];
coordlist = [xlist;ylist;zlist]';
xlist_Foci = [Results_micro.X];
ylist_Foci = [Results_micro.Y];
zlist_Foci = [Results_micro.Z];
coordlist_micro = [xlist_Foci;ylist_Foci;zlist_Foci]';
[k_micro,dist] = dsearchn(coordlist,coordlist_micro);
for i = 1:length(Results_all)
    fociindexlist = find(k_micro==i);
    Results_all(i).MN = fociindexlist;
end
%Map rupture MN
bafposmask = [Results_micro.BANF1] > 350;
indexlist = find(bafposmask);
for i = 1:length(Results_all)
    MNrupturemask(i) = 0;
    MNmask(i) = 0;
    if ~isempty(Results_all(i).MN)
        MNmask(i) = 1;
        for mni = 1:length(Results_all(i).MN)
            if any(indexlist==Results_all(i).MN)
                 MNrupturemask(i) = 1;        
            end
        end
    end
    Results_all(i).MNrupture = MNrupturemask(i);
    Results_all(i).MN = MNmask(i) ;
end
%Save Results
saveresult = struct2table(Results_all);
writetable(saveresult,[basedir '3D_all.csv'])

%% PN Rupture Comparison

mainmarkers = {'BANF1',...
           'NFkB','pH2AX','IRF3',...
            'SOX9',...
           'ELAVL4',...
          'SOX10',...
           'HOPX',...
           'Ki67',...
           'HIF1a','NDRG1',...
           'MX1'...
            };  
errlow = [];
diffperc = [];
sox2mask = [Results_all.SOX2mask];
PNrupturemask = [Results_all.FociN]>0;
for i = 1:length(mainmarkers)
    group = log10([Results_all.(mainmarkers{i})]./([Results_all.Hoechst1]).*mean([Results_all.Hoechst1])+1);% ;
    grouppos = group(PNrupturemask & sox2mask); %& sox2mask
    groupneg = group(~PNrupturemask & sox2mask);
    [~,p(i),est(i,1:2),~] = ttest2(grouppos,groupneg);
    diffperc(i) = mean(grouppos)/mean(groupneg)-1;
    error = est(i,1:2)/mean(groupneg)+1;
    errlow(i) = mean(grouppos)/mean(groupneg) - error(1);
end
cf = figure;
hold on
y = diffperc;
index = [];
[y,index] = sort(y);
xsorted = mainmarkers;
xsorted = xsorted(index);
x = categorical(xsorted);
x = reordercats(x,xsorted);
errlow = errlow(index);
p = p(index);
barh(x,y);
cf.Position = [1000 200 700 700];
er = errorbar(y,x,errlow,'.',"horizontal");
er.LineWidth = 1.5;
er.Color = 'Blue';
for i = 1:length(x)
    if p(i) < 0.01
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),"**",'Color','black','FontSize',12)
    elseif  p(i) < 0.05
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),'*','Color','black','FontSize',12)
    end
end
xlim([min(y)-1,max(y)+1])
text(0.1,x(2),"*p < 0.05","FontSize",12)
text(0.1,x(1),"**p < 0.01","FontSize",12)
ctitle = 'Difference in Primary Nucleus with Rupture versus no Rupture';
title(ctitle)
xlim([-0.2 0.2])
xlabel('Normalized Mean Nuclear Intensity Ratio')
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.pdf'],'Resolution',600)

%lamin A/C
f = figure;
f.Position = [500 500 300 200];
marker = 'LAMINAC';
y = [Results_morphology.(marker)];
grouppos = y([Results_morphology.Rupture] == 1);
grouppos = grouppos(grouppos>1);
groupneg = y([Results_morphology.Rupture] == 0);
groupneg = groupneg(groupneg>1);
[~,p] = ttest2(grouppos,groupneg)
boxdata = [groupneg grouppos];
boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
            repmat({'BAF+'},length(grouppos),1)];
%boxplot(boxdata,boxgroup)
hold on
    ksdensity(log10(grouppos+1))
    ksdensity(log10(groupneg+1))
    %legend('BAF+','BAF-')
hold off
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' marker '_PNrupture.pdf'],'Resolution',600)

%lamin B1
f = figure;
f.Position = [500 500 300 200];
marker = 'LaminB1';
y = [Results_morphology.(marker)];
grouppos = y([Results_morphology.Rupture] == 1);
grouppos = grouppos(grouppos>1);
groupneg = y([Results_morphology.Rupture] == 0);
groupneg = groupneg(groupneg>1);
[~,p] = ttest2(grouppos,groupneg)
boxdata = [groupneg grouppos];
boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
            repmat({'BAF+'},length(grouppos),1)];
%boxplot(boxdata,boxgroup)
xlim([2 4])
hold on
    ksdensity(log10(grouppos+1))
    ksdensity(log10(groupneg+1))
    %legend('BAF+','BAF-')
hold off
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' marker '_PNrupture.pdf'],'Resolution',600)

%DCX
f = figure;
f.Position = [500 500 300 200];
marker = 'DCX';
y = [Results_morphology.(marker)];
grouppos = y([Results_morphology.Rupture] == 1);
grouppos = grouppos(grouppos>1);
groupneg = y([Results_morphology.Rupture] == 0);
groupneg = groupneg(groupneg>1);
[~,p] = ttest2(grouppos,groupneg)
boxdata = [groupneg grouppos];
boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
            repmat({'BAF+'},length(grouppos),1)];
%boxplot(boxdata,boxgroup)
xlim([0 5])
hold on
    ksdensity(log10(grouppos+1))
    ksdensity(log10(groupneg+1))
    %legend('BAF+','BAF-')
hold off
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' marker '_PNrupture.pdf'],'Resolution',600)

%% Plot Generation
Metrics = {'Volume','Sphericity','Ellipticityoblate','Ellipticityprolate'};
sox2mask = [Results_morphology.SOX2mask]' == 1;
bafposmask = [Results_morphology.Rupture]' == 1;
for i = 1:4
    group = [Results_morphology.(Metrics{i})];
    grouppos = group(bafposmask & sox2mask);
    groupneg = group(~bafposmask & sox2mask);
    boxdata = [groupneg grouppos];
    boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
                repmat({'BAF+'},length(grouppos),1)];

   figure,
   boxplot(boxdata,boxgroup)
   %set(gca,'ytick',[])
   if i == 2
       %ylim([0 1])
   end
   %set(gca,'xtick',[])
    %violin({groupneg,grouppos},'xlabel',{'BAF-','BAF+'})
    %title(Metrics{i})
   %ttest
    [~,p,est,~] = ttest2(grouppos,groupneg);
    diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
    if p < 1
        if mean(est) > 0
            text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
                num2str(diffperc,2) '% Higher' newline ...
                'in BAF+'],'Color','black','FontSize',10)
        else
            text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
                num2str(diffperc,2) '% Lower' newline ...
                'in BAF+'],'Color','black','FontSize',10)
        end
    end
 
% %     figure,
% %     hold on
% %     ksdensity(groupneg)
% %     ksdensity(grouppos)
% %     hold off
% %     title(Metrics{i})
% %     ylabel('Probability')

    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' Metrics{i} '_Rupture.pdf'],'Resolution',600)
end
% 
% % Morphology Score
% IrregularIndex = [Results_morphology.IrregularIndex];
% sox2mask = [Results_morphology.SOX2mask]' == 1;
% bafposmask = [Results_morphology.Rupture]' == 1;
% 
% group = IrregularIndex;
% grouppos = group(bafposmask & sox2mask);
% groupneg = group(~bafposmask & sox2mask);
% boxdata = [groupneg grouppos];
% boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
%             repmat({'BAF+'},length(grouppos),1)];
% 
% figure,
% boxplot(boxdata,boxgroup)
% %set(gca,'ytick',[])
% %ylim([0 1])
% %set(gca,'xtick',[])
% %violin({groupneg,grouppos},'xlabel',{'BAF-','BAF+'})
% title('Morphology Irregularity Index')
% %ttest
% [~,p,est,~] = ttest2(grouppos,groupneg);
% diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
% if p < 1
%     if mean(est) > 0
%         text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
%             num2str(diffperc,2) '% Higher' newline ...
%             'in BAF+'],'Color','black','FontSize',10)
%     else
%         text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
%             num2str(diffperc,2) '% Lower' newline ...
%             'in BAF+'],'Color','black','FontSize',10)
%     end
% end
% exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/MII_Rupture_boxplot.jpg'],'Resolution',600)

Metrics = {'Volume','Sphericity','Ellipticityoblate','Ellipticityprolate'};
sox2mask = [Results_morphology.SOX2mask]' == 1;
bafposmask = [Results_morphology.Rupture]' == 1;
for i = 1:2
f = figure;
f.Position = [500 500 300 200];
hold on
group = [Results_morphology.(Metrics{i})];
posmask = group>0.5;
grouppos = group(bafposmask & sox2mask & posmask');
groupneg = group(~bafposmask & sox2mask & posmask');
ksdensity(groupneg,'Support','positive')
ksdensity(grouppos,'Support','positive')
legend('BAF-','BAF+')
[~,p,est,~] = ttest2(grouppos,groupneg)
hold off
title(Metrics{i})
ylabel('Probability')
%Figure 6 G
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Densityplot_Rupture_' Metrics{i} '.pdf'],'Resolution',600)
end

%% Lineage morphology 3D
Metrics = {'Volume','Sphericity','Ellipticityoblate','Ellipticityprolate','IrregularIndex'};
celltypes = {'AC','MES','OPC','NPC'};
celltypes_markers = {'SOX9','NDRG1','SOX10','ELAVL4'};
celltypes_markers_cutoff = [371,600,800,200];
sox2mask = [Results_morphology.SOX2mask]' == 1;
for typei=1:4
    celltypemean(typei) = mean([Results_morphology.(celltypes_markers{typei})]);
end
for cellid = 1:length(Results_morphology)
    for typei = 1:4
        celltypescore(typei) = Results_morphology(cellid).(celltypes_markers{typei}) ...
            ./ celltypes_markers_cutoff(typei);
%         if typei == 2
%             mean1 = Results_morphology(cellid).('NDRG1') ./ mean([Results_morphology.('NDRG1')]);
%             mean2 = Results_morphology(cellid).('HIF1a') ./ mean([Results_morphology.('HIF1a')]);
%             celltypescore(typei) = mean([mean1,mean2]);
%         end
    end
    [~,assigntypei] = max(celltypescore);
    celltypeslist(cellid) = assigntypei; 
%     for typei = 1:4
%         if Results_morphology(cellid).(celltypes_markers{typei}) > celltypes_markers_cutoff(typei)
%             celltypeslist(cellid) = typei;
%         end
%     end
    [~,~,ix] = unique(celltypeslist);
    C = accumarray(ix,1).';
end
for typei = 1:4
    %posmask = [Results.(celltypes_markers{typei})]'>celltypes_markers_cutoff(typei);
    posmask = celltypeslist' == typei;
    f = figure(typei);
    for i = 1:4
        group = [Results_morphology.(Metrics{i})];
        grouppos = group(posmask & sox2mask);
        groupneg = group(~posmask & sox2mask);
        boxdata = [groupneg grouppos];
        boxgroup = [repmat({['Non-' celltypes{typei}]},length(groupneg),1); ...
                    repmat({[celltypes{typei}]},length(grouppos),1)];
    
        subplot(2,2,i)
        boxplot(boxdata,boxgroup)
        title(Metrics{i})
    
        %ttest
    
        [~,p,est,~] = ttest2(grouppos,groupneg);
        diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
        if p < 1
            if mean(est) > 0
                text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
                    num2str(diffperc,2) '% Higher' newline ...
                    'in ' celltypes{typei}],'Color','black','FontSize',10)
            else
                text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
                    num2str(diffperc,2) '% Lower' newline ...
                    'in ' celltypes{typei}],'Color','black','FontSize',10)
            end
        end
    
        f.Position = [500 100 1000 800]; 
       
    end
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' celltypes{typei} '_Morphology.pdf'],'Resolution',600)
end
ylimlist = {[360 400],[0.8 0.85],[], [0.4 0.45]};
for i = [1,2,4]
    figure,
    cmetrics = Metrics{i};
    data = [];
    cat = [];
    for typei = 1:4
        cmask = celltypeslist' == typei;
        cdata = [Results_morphology.(cmetrics)];
        cdata =  cdata(cmask& sox2mask);
        data = [data cdata ];
        cat = [cat repmat({[celltypes{typei}]},1,length(cdata))];
        catcell{typei} = cdata;
        barmean(typei) =  mean(cdata);
    end
    boxplot(data,cat)
   %vplot = violin(catcell,'xlabel',celltypes);
   % bar(categorical(celltypes),barmean)
    title(Metrics{i})
    ylim(ylimlist{i})
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' cmetrics '_Lineage.pdf'],'Resolution',600)
end

for i = [1,2,4]
    
    f = figure;
    %f.Position = [1000 500 700 400];
    hold on
    cmetrics = Metrics{i};
    data = [];
    cat = [];
    for typei = 1:4
        cmask = celltypeslist' == typei;
        cdata = [Results_morphology.(cmetrics)];
        posmask = cdata'>0;
        if i == 5
            posmask = cdata'>0 & cdata'<1;
        end
        cdata_pos =  cdata(cmask& sox2mask & posmask);
        data = [data cdata_pos ];
        violindata{typei} = cdata_pos; 
        cat = [cat repmat({[celltypes{typei}]},1,length(cdata_pos))];

        grouppos = cdata_pos;
        groupneg = cdata(~cmask& sox2mask & posmask);
        [~,p,est,~] = ttest2(grouppos,groupneg);
        diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
        if p < 0.1
            if mean(est) > 0
                text(typei+0.1,mean(grouppos)*1.2,['p = ' num2str(p) newline ... 
                    num2str(diffperc,2) '% Higher' newline ...
                    'in ' celltypes{typei}],'Color','black','FontSize',9)
            else
                text(typei+0.1,mean(grouppos)*1.2,['p = ' num2str(p) newline ...
                    num2str(diffperc,2) '% Lower' newline ...
                    'in ' celltypes{typei}],'Color','black','FontSize',9)
            end
        end
    end
    boxplot(data,cat)
    xlabel = celltypes;
    %violin(violindata,'xlabel',xlabel)
    title(cmetrics)
    hold off
    
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Lineage_all_' cmetrics '.pdf'],'Resolution',600)
end
%% MN VS PN



Metrics = {'Volume','Sphericity','Ellipticityoblate','Ellipticityprolate'};
% Metrics = {'TREX1','pH2AX','cGAS','Ki67','DAPI'};

for i = 1:4
    grouppos = [Results_micro.(Metrics{i})];
    grouppos = grouppos([Results_micro.SOX2]>0);
    groupneg = [Results_morphology.(Metrics{i})];
    groupneg = groupneg([Results_morphology.SOX2]>300);
    boxdata = [groupneg grouppos];
    boxgroup = [repmat({'Primary'},length(groupneg),1); ...
                repmat({'Micro'},length(grouppos),1)];

     figure,
    boxplot(boxdata,boxgroup)
%     title(['Micronucleus ' Metrics{i} ' Ruptured vs Unruptured' ])
%     ylabel(Metrics{i})
%ylim([0.5 1])
%set(gca,'ytick',[])
    %ttest
%     Vdata{1} = groupneg;
%     Vdata{2} = grouppos;
%     violin(Vdata,'xlabel',{'BAF-','BAF+'})
%     scatter([grouppos,repmat(1,4)])
%     scatter([groupneg,1)
%     title(['Micronucleus ' Metrics{i} ' Ruptured vs Unruptured' ])
%     ylabel(Metrics{i})
% 
    [~,p,est,~] = ttest2(grouppos,groupneg);
    diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
    if p < 1
        if mean(est) > 0
            text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
                num2str(diffperc) '% Higher' newline ...
                'in BAF+'],'Color','black','FontSize',10)
        else
            text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
                num2str(diffperc) '% Lower' newline ...
                'in BAF+'],'Color','black','FontSize',10)
        end
    end
end

%% Micronucelus Morphology

bafposmask = [Results_micro.BANF1] > 300;
sox2mask = [Results_micro.SOX2] > 0;

Metrics = {'Volume','Sphericity','Ellipticityoblate','Ellipticityprolate'};
% Metrics = {'TREX1','pH2AX','cGAS','Ki67','DAPI'};

for i = 1:4
    group = [Results_micro.(Metrics{i})];
    grouppos = group(bafposmask & sox2mask);
    groupneg = group(~bafposmask & sox2mask);
    boxdata = [groupneg grouppos];
    boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
                repmat({'BAF+'},length(grouppos),1)];

     figure,
    boxplot(boxdata,boxgroup)
%     title(['Micronucleus ' Metrics{i} ' Ruptured vs Unruptured' ])
%     ylabel(Metrics{i})
if i == 2
%ylim([0 1])
end
%set(gca,'ytick',[])
    %ttest
%     Vdata{1} = groupneg;
%     Vdata{2} = grouppos;
%     violin(Vdata,'xlabel',{'BAF-','BAF+'})
%     scatter([grouppos,repmat(1,4)])
%     scatter([groupneg,1)
%     title(['Micronucleus ' Metrics{i} ' Ruptured vs Unruptured' ])
%     ylabel(Metrics{i})
% 
    [~,p,est,~] = ttest2(grouppos,groupneg);
    diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
    if p < 1
        if mean(est) > 0
            text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
                num2str(diffperc) '% Higher' newline ...
                'in BAF+'],'Color','black','FontSize',10)
        else
            text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
                num2str(diffperc) '% Lower' newline ...
                'in BAF+'],'Color','black','FontSize',10)
        end
    end

    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' Metrics{i} '_MN_Rupture.pdf'],'Resolution',600)
end

bafposmask = [Results_micro.BANF1] > 300;
sox2mask = [Results_micro.SOX2] > 0;
for i = 1:2
f = figure;
f.Position = [500 500 300 200];
hold on
group = [Results_micro.(Metrics{i})];
posmask = group>0.5;
grouppos = group(bafposmask & sox2mask & posmask);
groupneg = group(~bafposmask & sox2mask & posmask);
ksdensity(groupneg,'Support','positive')
ksdensity(grouppos,'Support','positive')
legend('BAF-','BAF+')
if i == 1
    xlim([-5 20]);
elseif i == 2
    xlim([0.6 1.1])
end
hold off
title(Metrics{i})
ylabel('Probability')
[~,p,est,~] = ttest2(grouppos,groupneg)
%Figure 6 G
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Densityplot_Rupture_MN_' Metrics{i} '.pdf'],'Resolution',600)
end

bafposmask = [Results_micro.BANF1] > 300;
sox2mask = [Results_micro.SOX2] > 0;
for i = 1:2
f = figure;
f.Position = [500 500 300 200];
hold on
group = [Results_micro.(Metrics{i})];
posmask = group>0.5;
grouppos = group(bafposmask & sox2mask & posmask);
histogram(grouppos)
%legend('BAF-','BAF+')
if i == 1
    xlim([0 25]);
elseif i == 2
    xlim([0.7 1])
end
hold off
title(Metrics{i})
ylabel('Count')
[~,p,est,~] = ttest2(grouppos,groupneg)
%Figure 6 G
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Densityplot_MN_' Metrics{i} '.pdf'],'Resolution',600)
end
%% Associate MN

sox2mask = [Results_all.SOX2]' > 200;
%Metrics = {'Volume','Sphericity','Ellipticityoblate','Ellipticityprolate'};
Metrics = {'TREX1','pH2AX','cGAS','Ki67','Hoechst1'};
% 
% for i = 1:5
%     group = log10([Results_all.(Metrics{i})]+1.1)-log10([Results_all.Hoechst1]+1.1);
%     outliermask = (group < 0.5)';
%     grouppos = group(MNmask' & sox2mask & outliermask); %& sox2mask
%     groupneg = group(~MNmask' & sox2mask & outliermask);
%     boxdata = [groupneg grouppos];
%     boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
%                 repmat({'BAF+'},length(grouppos),1)];
% 
%     %ttest
%     figure,
%     Vdata{1} = groupneg;
%     Vdata{2} = grouppos;
%     violin(Vdata,'xlabel',{'Without MN','With MN'})
% %     scatter([grouppos,repmat(1,4)])
% %     scatter([groupneg,1)
%     title(['Primary Nucleus with and without MN' ])
%     ylabel(Metrics{i})
% 
%     [~,p,est,~] = ttest2(grouppos,groupneg);
%     diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
%     if p < 1
%         if mean(est) > 0
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
%                 num2str(diffperc) '% Higher' newline ...
%                 'in BAF+'],'Color','black','FontSize',10)
%         else
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
%                 num2str(diffperc) '% Lower' newline ...
%                 'in BAF+'],'Color','black','FontSize',10)
%         end
%     end
%    % set(gca, 'YScale', 'log')
%     exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' Metrics{i} '_PNwithMN.jpg'],'Resolution',600)
% end
% 
% for i = 1:5
%     group = log10([Results_all.(Interferons{i})]+1.1)-log10([Results_all.Hoechst1]+1.1);
%     grouppos = group(MNmask' & sox2mask); %& sox2mask
%     groupneg = group(~MNmask' & sox2mask);
%     boxdata = [groupneg grouppos];
%     boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
%                 repmat({'BAF+'},length(grouppos),1)];
% 
%     %ttest
%     figure,
%     Vdata{1} = groupneg;
%     Vdata{2} = grouppos;
%     violin(Vdata,'xlabel',{'Without MN','With MN'})
% %     scatter([grouppos,repmat(1,4)])
% %     scatter([groupneg,1)
%     title(['Primary Nucleus with and without MN' ])
%     ylabel(Interferons{i})
% 
%     [~,p,est,~] = ttest2(grouppos,groupneg);
%     diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
%     if p < 1
%         if mean(est) > 0
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
%                 num2str(diffperc) '% Higher' newline ...
%                 'in BAF+'],'Color','black','FontSize',10)
%         else
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
%                 num2str(diffperc) '% Lower' newline ...
%                 'in BAF+'],'Color','black','FontSize',10)
%         end
%     end
%    % set(gca, 'YScale', 'log')
%     exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' Interferons{i} '_PNwithMN.jpg'],'Resolution',600)
% end

mainmarkers = {'BANF1',...
           'NFkB','pH2AX','IRF3',...
            'SOX9',...
           'ELAVL4',...
          'SOX10',...
           'HOPX',...
           'Ki67',...
           'HIF1a','NDRG1',...
           'MX1'...
            }; 
p = [];
 diffperc = [];
for i = 1:length(mainmarkers)
    group = ([Results_all.(mainmarkers{i})]+1.1) ./ ([Results_all.Hoechst1]+1.1);
    grouppos = group(MNmask' & sox2mask); %& sox2mask
    groupneg = group(~MNmask' & sox2mask);
    [~,p(i),est(i,1:2),~] = ttest2(grouppos,groupneg);
    diffperc(i) = mean(grouppos)/mean(groupneg)-1;
    error = est(i,1:2)/mean(groupneg)+1;
    errlow(i) = mean(grouppos)/mean(groupneg) - error(1);
end
cf = figure;
hold on
y = diffperc;
index = [];
[y,index] = sort(y);
xsorted = mainmarkers;
xsorted = xsorted(index);
x = categorical(xsorted);
x = reordercats(x,xsorted);
errlow = errlow(index);
p = p(index);
barh(x,y);
cf.Position = [1000 200 540 700];
er = errorbar(y,x,errlow,'.',"horizontal");
er.LineWidth = 1.5;
er.Color = 'Blue';
for i = 1:length(x)
    if p(i) < 0.01
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),"**",'Color','black','FontSize',12)
    elseif  p(i) < 0.05
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),'*','Color','black','FontSize',12)
    end
end
xlim([min(y)-0.3,max(y)+0.3])
text(0.4,x(2),"*p < 0.05","FontSize",12)
text(0.4,x(1),"**p < 0.01","FontSize",12)
ctitle = 'Difference in Primary Nucleus with Micronucleus versus no Micronucleus';
title(ctitle)
xlabel('Normalized Mean Nuclear Intensity Ratio')
exportgraphics(cf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.pdf'],'Resolution',600)
%% Associate microrupture

sox2mask = [Results_all.SOX2]' > 200;
%Metrics = {'Volume','Sphericity','Ellipticityoblate','Ellipticityprolate'};
Metrics = {'TREX1','pH2AX','cGAS','Ki67','DAPI'};
% 
% for i = 1:4
%     group = log10([Results_all.(Metrics{i})]./([Results_all.Hoechst1]).*mean([Results_all.Hoechst1])+1);
%     grouppos = group(MNrupturemask' & MNmask' & sox2mask); %& sox2mask
%     groupneg = group(~MNrupturemask' & sox2mask);
%     boxdata = [groupneg grouppos];
%     boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
%                 repmat({'BAF+'},length(grouppos),1)];
% 
%      figure,
%     %ttest
%     Vdata{1} = groupneg;
%     Vdata{2} = grouppos;
%     violin(Vdata,'xlabel',{'BAF- MN' 'BAF+ MN'})
% %     scatter([grouppos,repmat(1,4)])
% %     scatter([groupneg,1)
%     title(['Primary Nucleus with vs without MN BAF+ ' ])
%     ylabel(Metrics{i})
% 
%     [~,p,est,~] = ttest2(grouppos,groupneg);
%     diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
%     if p < 1
%         if mean(est) > 0
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
%                 num2str(diffperc) '% Higher' newline ...
%                 'in BAF+'],'Color','black','FontSize',10)
%         else
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
%                 num2str(diffperc) '% Lower' newline ...
%                 'in BAF+'],'Color','black','FontSize',10)
%         end
%     end
%    % set(gca, 'YScale', 'log')
%     exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' Metrics{i} '_PNwithMNRupture.jpg'],'Resolution',600)
% end
% 
% for i = 1:4
%     group = log10([Results_all.(Interferons{i})]./([Results_all.Hoechst1]).*mean([Results_all.Hoechst1])+1);
%     grouppos = group(MNrupturemask' & MNmask' & sox2mask); %& sox2mask
%     groupneg = group(~MNrupturemask' & sox2mask);
%     boxdata = [groupneg grouppos];
%     boxgroup = [repmat({'BAF-'},length(groupneg),1); ...
%                 repmat({'BAF+'},length(grouppos),1)];
% 
%      figure,
%     %ttest
%     Vdata{1} = groupneg;
%     Vdata{2} = grouppos;
%     violin(Vdata,'xlabel',{'BAF- MN' 'BAF+ MN'})
%     title(['Primary Nucleus with vs without MN BAF+ ' ])
%     ylabel(Interferons{i})
% 
%     [~,p,est,~] = ttest2(grouppos,groupneg);
%     diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
%     if p < 1
%         if mean(est) > 0
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
%                 num2str(diffperc) '% Higher' newline ...
%                 'in BAF+'],'Color','black','FontSize',10)
%         else
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
%                 num2str(diffperc) '% Lower' newline ...
%                 'in BAF+'],'Color','black','FontSize',10)
%         end
%     end
%    % set(gca, 'YScale', 'log')
%     exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' Interferons{i} '_PNwithMNRupture.jpg'],'Resolution',600)
% end

mainmarkers = {...
           'NFkB','pH2AX','IRF3',...
            'SOX9',...
           'ELAVL4',...
          'SOX10',...
           'HOPX',...
           'Ki67',...
           'HIF1a','NDRG1',...
           'MX1'...
            };  
errlow = [];
diffperc = [];
for i = 1:length(mainmarkers)
    group = log10([Results_all.(mainmarkers{i})]./([Results_all.Hoechst1]).*mean([Results_all.Hoechst1])+1);% ;
    %group = [Results_all.(mainmarkers{i})]./([Results_all.Hoechst1]).*mean([Results_all.Hoechst1]);
    grouppos = group(MNrupturemask' & MNmask' & sox2mask); %& sox2mask
    groupneg = group(~MNrupturemask' & sox2mask);
    [~,p(i),est(i,1:2),~] = ttest2(grouppos,groupneg);
    diffperc(i) = mean(grouppos)/mean(groupneg)-1;
    error = est(i,1:2)/mean(groupneg)+1;
    errlow(i) = mean(grouppos)/mean(groupneg) - error(1);
end
cf = figure;
hold on
y = diffperc;
index = [];
[y,index] = sort(y);
xsorted = mainmarkers;
xsorted = xsorted(index);
x = categorical(xsorted);
x = reordercats(x,xsorted);
errlow = errlow(index);
p = p(index);
barh(x,y);
cf.Position = [1000 200 700 700];
er = errorbar(y,x,errlow,'.',"horizontal");
er.LineWidth = 1.5;
er.Color = 'Blue';
for i = 1:length(x)
    if p(i) < 0.01
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),"**",'Color','black','FontSize',12)
    elseif  p(i) < 0.05
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),'*','Color','black','FontSize',12)
    end
end
xlim([min(y)-1,max(y)+1])
text(1,x(2),"*p < 0.05","FontSize",12)
text(1,x(1),"**p < 0.01","FontSize",12)
ctitle = 'Difference in Primary Nucleus with Micronucleus Rupture versus no Micronucleus';
title(ctitle)
xlim([-1 2.3])
xlabel('Normalized Mean Nuclear Intensity Ratio');
exportgraphics(cf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.pdf'],'Resolution',600)

mainmarkers = {...
           'NFkB','pH2AX','IRF3',...
            'SOX9',...
           'ELAVL4',...
          'SOX10',...
           'HOPX',...
           'Ki67',...
           'HIF1a','NDRG1',...
           'MX1'...
            };  

% mainmarkers = {'BANF1','STING','TREX1'...
%            'NFkB','pH2AX','IRF3'...
%            'SOX2','cGAS','SOX9'...
%            'ELAVL4','DCX','OLIG2'...
%            'EGFR','CEBPB','SOX10'...
%            'TMEM119','LAMINAC','HOPX'...
%            'Ki67','LaminB1','LaminB2'...
%            'HIF1a','NDRG1','IFITM1'...
%            'CD31','MX1','INSM1'...
%            'BANF1_2','cGAS_2','pH2AX_2'...
%             }; 
errlow = [];
p = [];
est = [];
diffperc = [];
for i = 1:length(mainmarkers)
    %group = log10([Results_all.(mainmarkers{i})]./([Results_all.Hoechst1]).*mean([Results_all.Hoechst1])+1);% ;
    group = [Results_all.(mainmarkers{i})]./([Results_all.Hoechst1]).*mean([Results_all.Hoechst1]);
    grouppos = group(MNrupturemask' & MNmask' & sox2mask); %
    groupneg = group(~MNrupturemask' & MNmask' & sox2mask);
    [~,p(i),est(i,1:2),~] = ttest2(grouppos,groupneg);
    diffperc(i) = mean(grouppos)/mean(groupneg)-1;
    error = est(i,1:2)/mean(groupneg)+1;
    errlow(i) = mean(grouppos)/mean(groupneg) - error(1);
end
cf = figure;
hold on
y = diffperc;
index = [];
[y,index] = sort(y);
xsorted = mainmarkers;
xsorted = xsorted(index);
x = categorical(xsorted);
x = reordercats(x,xsorted);
errlow = errlow(index);
p = p(index);
barh(x,y);
cf.Position = [1000 200 700 700];
er = errorbar(y,x,errlow,'.',"horizontal");
er.LineWidth = 1.5;
er.Color = 'Blue';
for i = 1:length(x)
    if p(i) < 0.01
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),"**",'Color','black','FontSize',12)
    elseif  p(i) < 0.05
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),'*','Color','black','FontSize',12)
    end
end
xlim([min(y)-1,max(y)+1])
text(1,x(2),"*p < 0.05","FontSize",12)
text(1,x(1),"**p < 0.01","FontSize",12)
ctitle = 'Difference in Primary Nucleus with Micronucleus Rupture versus no Micronucleus Rupture';
title(ctitle)
xlim([-1 1.8])
xlabel('Normalized Mean Nuclear Intensity Ratio')
exportgraphics(cf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.pdf'],'Resolution',600)

%% NPC and OPC specific graph
sox2mask = [Results_all.SOX2]' > 200;
typei = 4;
i = 1;
posmask = [Results.(celltypes_markers{typei})]'>celltypes_markers_cutoff(typei);
group = [Results.(Metrics{i})];
grouppos = group(posmask);
groupneg = group(~posmask);
f = figure(1);
hold on
ksdensity(grouppos)
ksdensity(groupneg)
legend(celltypes{typei},['Non-' celltypes{typei}])
hold off
[~,p,~,~] = ttest2(grouppos,groupneg);
title(['p = ' num2str(p) ' ' Metrics{i} ' '...
                    num2str((mean(grouppos)-mean(groupneg))./mean(groupneg)*100,2) ...
                    '% Higher ' ...
                    'in ' celltypes{typei} ' [' celltypes_markers{typei} ']'] ...
                    ,'Color','black','FontSize',10)
xlabel('Volume (um^3)')
ylabel('Probability')

exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' celltypes{typei} '_' Metrics{i} '.jpg'],'Resolution',600)
typei = 3;
i = 2;
posmask = [Results.(celltypes_markers{typei})]'>celltypes_markers_cutoff(typei);%
group = [Results.(Metrics{i})];
grouppos = group(posmask&sox2mask);
groupneg = group(~posmask&sox2mask);
f = figure(2);
hold on
ksdensity(grouppos)
ksdensity(groupneg)
legend(celltypes{typei},['Non-' celltypes{typei}])
hold off
[~,p,~,~] = ttest2(grouppos,groupneg);
title(['p = ' num2str(p) ' ' Metrics{i} ' '...
                    num2str((mean(grouppos)-mean(groupneg))./mean(groupneg)*100,2) ...
                    '% Higher ' ...
                    'in ' celltypes{typei} ' [' celltypes_markers{typei} ']'] ...
                    ,'Color','black','FontSize',10)
xlabel('Volume (um^3)')
ylabel('Probability')
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' celltypes{typei} '_' Metrics{i} '.jpg'],'Resolution',600)

%% Frequency of Rupture
sox2mask = [Results_all.SOX2]' > 200;
profile on
load("C:\Users\brc642\Desktop\temp.mat")
data_2d = data;
data = [Results_all.FociN];
data = data(sox2mask);
lambdahat = poissfit(data);
x_pois = 0:7;
y_pois = poisspdf(x_pois,lambdahat);
f = figure;
    f.Position = [400 400 500 500];
    hold on
    histogram1 = histogram(data,'BinMethod','integers');
%     histogram1.Normalization = 'probability';
%     histogram2 = histogram(data_2d,'BinMethod','integers');
%     histogram2.Normalization = 'probability';
    legend('Frequency Distribution 3D','Frequency Distribution 2D')
    xlim([-0.5 7.5])
    ylim([0 4000])
    title('Rupture Frequency')
    xlabel('Number of Rupture per Primary Nucleus')
     ylabel('Percentage')
    plot1 = plot(x_pois,y_pois*sum(sox2mask));
    legend('Frequency Distribution',['Poisson Estimation (λ = ' num2str(lambdahat) ')'])
    hold off
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/BAF_Rupturerate.pdf'],'Resolution',600)

%Markers
Channel = 'LAMINAC';
sox2mask = [Results_all.SOX2]' > 200;
ycells = [];
for xi = 0:5
    dataraw= [Results_all(sox2mask).(Channel)];
    data = log10(dataraw);
    mask = [Results_all(sox2mask).FociN] == xi;
    posmask = data>0;
    ycells{xi+1} = data(mask&posmask);
end 
figure,
    xlabellist = {'0','1','2','3','4','5'};
    violin(ycells,'xlabel',xlabellist)
    xlabel('Number of Rupture');
    ylabel(['log10(' Channel ')'])
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Lamin_RuptureRate.jpg'],'Resolution',600)
Channel = 'LAMINAC';
sox2mask = [Results_all.SOX2]' > 200;
ycells = [];
for xi = 0
    dataraw= [Results_all(sox2mask).(Channel)];
    data = log10(dataraw);
    mask = [Results_all(sox2mask).FociN] == xi;
    posmask = data>0;
    ycells{xi+1} = data(mask&posmask);
end 
for xi = 1
    dataraw= [Results_all(sox2mask).LAMINAC];
    data = log10(dataraw);
    mask = [Results_all(sox2mask).FociN] >= xi;
    posmask = data>0;
    ycells{xi+1} = data(mask&posmask);
    mean(ycells{xi+1});
end
figure,
    xlabellist = {'0','1+'};
    violin(ycells,'xlabel',xlabellist)
    xlabel('Number of Rupture');
    ylabel(['log10(' Channel ')'])
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Lamin_RuptureRate.jpg'],'Resolution',600)

%Markers
Channel = 'Volume';
sox2mask = [Results_all.SOX2]' > 200;
ycells = [];
for xi = 0:3
    dataraw= [Results_all(sox2mask).(Channel)];
    data = log10(dataraw);
    mask = [Results_all(sox2mask).FociN] == xi;
    posmask = data>0&data<prctile(data,95)&data>prctile(data,5);
    ycells{xi+1} = data(mask&posmask);
end 
for xi = 4
    dataraw= [Results_all(sox2mask).(Channel)];
    data = log10(dataraw);
    mask = [Results_all(sox2mask).FociN] >= xi;
    posmask = data>0&data<prctile(data,95)&data>prctile(data,5);
    ycells{xi+1} = data(mask&posmask);
end
f = figure;
f.Position = [500 500 500 500]
    xlabellist = {'0','1','2','3','4+'};
    violin(ycells,'xlabel',xlabellist)
    xlabel('Number of Rupture');
    ylabel(['log10(' Channel ')'])
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Volume_RuptureRate.pdf'],'Resolution',600)

%Markers
Channel = 'Sphericity';
sox2mask = [Results_all.SOX2]' > 200;
ycells = [];
for xi = 0:3
    dataraw= [Results_all(sox2mask).(Channel)];
    data = dataraw;
    mask = [Results_all(sox2mask).FociN] == xi;
    posmask = data>0&data<prctile(data,95)&data>prctile(data,5);
    ycells{xi+1} = data(mask&posmask);
end 
for xi = 4
    dataraw= [Results_all(sox2mask).(Channel)];
    data = dataraw;
    mask = [Results_all(sox2mask).FociN] >= xi;
    posmask = data>0&data<prctile(data,95)&data>prctile(data,5);
    ycells{xi+1} = data(mask&posmask);
end
f = figure;
f.Position = [500 500 500 500]
    xlabellist = {'0','1','2','3','4+'};
    violin(ycells,'xlabel',xlabellist)
    xlabel('Number of Rupture');
    ylabel([' Channel '])
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Sphericity_RuptureRate.pdf'],'Resolution',600)


%Markers
Channel = 'Volume';
sox2mask = [Results_all.SOX2]' > 200;
ycells = [];
for xi = 0:4
    dataraw= [Results_all(sox2mask).(Channel)];
    data = log10(dataraw);
    mask = [Results_all(sox2mask).FociN] == xi;
    posmask = data>0;
    ycells{xi+1} = data(mask&posmask);
end 
figure,
    xlabellist = {'0','1','2','3','4+'};
    violin(ycells,'xlabel',xlabellist)
    xlabel('Number of Rupture');
    ylabel(['log10(' Channel ')'])
    title('Rupture Frequency vs. Volume')
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Volume_RuptureRate.pdf'],'Resolution',600)

%Markers
Channel = 'pH2AX';
channelmask = [Results_all.(Channel)]' > 200;
ycells = [];
for xi = 0:5
    dataraw= [Results_all(channelmask).(Channel)];
    data = log10(dataraw);
    mask = [Results_all(channelmask).FociN] == xi;
    posmask = data>0;
    ycells{xi+1} = data(mask&posmask);
end 
figure,
    xlabellist = {'0','1','2','3','4','5'};
    violin(ycells,'xlabel',xlabellist)
    xlabel('Number of Rupture');
    ylabel(['log10(' Channel ')'])
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Volume_RuptureRate.jpg'],'Resolution',600)

%Markers
Channel = 'IrregularIndex';
channelmask = [Results_all.(Channel)]' > 0;
ycells = [];
for xi = 0:4
    dataraw= [Results_all(channelmask).(Channel)];
    %data = log10(dataraw);
    mask = [Results_all(channelmask).FociN] == xi;
    posmask = dataraw>0&dataraw<1;
    ycells{xi+1} = dataraw(mask&posmask);
    if xi>0
        [~,p] = ttest2(ycells{xi+1},ycells{xi})
    end
end 
    f=figure;
    f.Position = [1000 100 500 700];
    xlabellist = {'0','1','2','3','4+'};
    violin(ycells,'xlabel',xlabellist)
    xlabel('Number of Rupture');
    ylabel(['IrregularIndex'])
    title('Irregularity Index v. Frequency of Rupture in 3D GBM')
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/NRUP_IrregularIndex.pdf'],'Resolution',600)



%% Rupture Shape

Metrics = {'Volume','Sphericity','Ellipticityoblate','Ellipticityprolate'};
xdata = [Results_Foci.Volume];
xdata = log10(xdata+1);
hist = histogram(xdata);
hist.Normalization = 'probability';
xlabel('Rupture Size (um3)')
ylabel('Probability')
titlename = 'BAF Foci Volume';
title('BAF Foci Volume')
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' titlename '.pdf'],'Resolution',600)

metric = 'Sphericity';
xdata = [Results_Foci.(metric)];
hist = histogram(xdata);
hist.Normalization = 'probability';
xlabel(metric)
ylabel('Probability')
titlename = 'BAF Foci Sphericity';
title(['BAF Foci ' metric])
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' titlename '.pdf'],'Resolution',600)

hold on
xdata = [Results_Foci.pH2AX];
xdata = log10(xdata+1);
hist = histogram(xdata,'BinWidth',0.25);
hist.Normalization = 'probability';
xdata = [Results_all.pH2AX];
xdata = log10(xdata+1);
hist2 = histogram(xdata,'BinWidth',0.25);
hist2.Normalization = 'probability';
xlabel('log(pH2AX mean intensity level)')
ylabel('Probability')
titlename = 'BAF Foci pH2AX';
title(titlename)
legend('Foci','Nuclear')
hold off
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' titlename '.jpg'],'Resolution',600)

hold on
xdata = [Results_Foci.BANF1];
xdata = log10(xdata+1);
hist = histogram(xdata,'BinWidth',0.25);
hist.Normalization = 'probability';
xdata = [Results_all.pH2AX];
xdata = log10(xdata+1);
hist2 = histogram(xdata,'BinWidth',0.25);
hist2.Normalization = 'probability';
xlabel('log(BANF1 mean intensity level)')
ylabel('Probability')
titlename = 'BAF Foci BANF1';
title(titlename)
legend('Foci','Nuclear')
hold off
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' titlename '.jpg'],'Resolution',600)


xdata = [Results_Foci.pH2AX];
Fociposmask = xdata > 1000;
Focinegmask = ~Fociposmask;
data = [Results_Foci.Volume];
grouppos = data(Fociposmask);
grouppos = log10(grouppos+1);
groupneg = data(Focinegmask);
groupneg = log10(groupneg+1);
figure
    hold on
    violin({groupneg,grouppos},'xlabel',{'BAF-','BAF+'})
    hold off
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' titlename '.jpg'],'Resolution',600)

%% Rupture shape classification
Volumelist = [Results_Foci.Volume];
Sphericitylist = [Results_Foci.Sphericity];
Elliplicitylist = [Results_Foci.Ellipticityprolate];
range = prctile(Volumelist,[5 95]);
mask1 = Volumelist>range(1) & Volumelist<range(2); 
range = prctile(Sphericitylist,[5 95]);
mask2 = Sphericitylist>range(1) & Sphericitylist<range(2); 
range = prctile(Elliplicitylist,[5 95]);
mask3 = Elliplicitylist>range(1) & Elliplicitylist<range(2); 
mask123 = mask1 & mask2 & mask3;
n = [];
n(1,:) = (Volumelist- prctile(Volumelist,5)) ./ prctile(Volumelist,95);
n(2,:) = Sphericitylist;%(Sphericitylist- prctile(Sphericitylist,5))./ prctile(Sphericitylist,95);
n(3,:) = Elliplicitylist;%(Elliplicitylist- prctile(Elliplicitylist,5))./ prctile(Elliplicitylist,95);
nscore = n(1,:) ./ n(2,:);% .* n(3,:);
soxmask = [Results_Foci.SOX2] > 200;
% soxmask = soxmask&mask123;
%lowsphericity = Sphericitylist < 0.85;
irregularIndex = nscore;%Volumelist;%.*((1-Sphericitylist));
irregularFoci = irregularIndex>prctile(irregularIndex,75);
regularFoci = irregularIndex<prctile(irregularIndex,25); 
% irregularFoci = (Volumelist > 10) & (Sphericitylist < 0.9);
% figure,
%     hold on
%     hist = histogram(log10(Volumelist),'BinWidth',0.1);
%     hist.Normalization = "count";
%     hist2 = histogram(log10(Volumelist(irregularFoci)),'BinWidth',0.1);
%     hist2.Normalization = "count";
%     hold off
%     xlabel('log10(Rupture Size) (um3)')
%     ylabel('Count')
%     ctitle = 'BAF Foci Volume';
%     title(ctitle)
%     %set(gca, 'YScale', 'log')
%     legend('Normal Rupture','Irregular Rupture')
%     exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.jpg'],'Resolution',600)
% figure,
%     hold on
%     hist = histogram(Sphericitylist,'BinWidth',0.02);
%     hist.Normalization = "count";
%     hist2 = histogram(Sphericitylist(irregularFoci),'BinWidth',0.02);
%     hist2.Normalization = "count";
%     hold off
%     xlabel('Sphericity')
%     ylabel('Count')
%     ctitle = 'BAF Foci Sphericity';
%     title(ctitle)
%     legend('Normal Rupture','Irregular Rupture')
%     %set(gca, 'YScale', 'log')
%     exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.jpg'],'Resolution',600)figure,
figure,
    hold on
    hist = histogram(irregularIndex,'BinWidth',0.25);
    hist.Normalization = "count";
    hold off
    %xlabel('Morphology Irregular Index')
    ylabel('Count')
    ctitle = 'Morphology Irregularity Index';
    title(ctitle)
    xlim([-1 6])
    set(gca, 'YScale', 'log')
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/MorphologyIrregularIndex.pdf'],'Resolution',600)
figure,
    hold on
    hist = histogram(n(2,:));
    hist.Normalization = "count";
    hold off
    %xlabel('Morphology Irregular Index')
    ylabel('Count')
    ctitle = 'Normalized Sphericity Score';
    title(ctitle)
    set(gca, 'YScale', 'log')
    %exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/MorphologyIrregularIndex.pdf'],'Resolution',600)
est = [];
diffperc = [];
Metrics = {'pH2AX',Interferons{:}}; %'BANF1','LAMINAC','LaminB1','LaminB2',,'cGAS','Hoechst1'
for i = 1:length(Metrics)
    group = log10(([Results_Foci.(Metrics{i})])+1);%
    posmask = group > 0;
    grouppos = group(irregularFoci & soxmask & posmask);
    groupneg = group(regularFoci & soxmask & posmask);
    boxdata = [groupneg grouppos];

    f = figure;
    f.Position = [1000 500 200 400];
    violin({groupneg,grouppos},'xlabel',{'Regular','Irregular'})
    title(Metrics{i})
    ylabel(['Log10(' Metrics{i}  ')'])
%   set(gca, 'YScale', 'log')
   %ttest
    [~,p(i),est(i,1:2),~] = ttest2(grouppos,groupneg);
    diffperc(i) = mean(grouppos)-mean(groupneg);
    %error = est(i,1:2)/mean(groupneg)+1;
    %errlow(i) = mean(grouppos)/mean(groupneg) - error(1);
    errlow(i) = est(i,1);
%     if p < 1
%         if mean(est) > 0
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ... 
%                 num2str(diffperc,2) '% Higher' newline ...
%                 'in Irregular Foci'],'Color','black','FontSize',10)
%         else
%             text(1.3,mean(grouppos),['p = ' num2str(p) newline ...
%                 num2str(diffperc,2) '% Lower' newline ...
%                 'in Irregular Foci'],'Color','black','FontSize',10)
%         end
%     end
    titlename = [Metrics{i} 'Rupture Morphology'];
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/BAF_Rupture_Shape_' titlename '.pdf'],'Resolution',600)
end

% Metrics = {'BANF1','LAMINAC','LaminB1','LaminB2','pH2AX','cGAS','Hoechst1'};
% for i = 1:length(Metrics)
%     group = [Results_Foci.(Metrics{i})];
%     cmatrix = corr([group(soxmask)' irregularIndex(soxmask)']);
%     table{i,2} = cmatrix(2);
%     table{i,1} = Metrics{i};
% end

cf = figure;
hold on
y = diffperc;
index = [];
[y,index] = sort(y);
xsorted = Metrics;
xsorted = xsorted(index);
x = categorical(xsorted);
x = reordercats(x,xsorted);
errlow = errlow(index);
p = p(index);
barh(x,y);
cf.Position = [1000 200 300 300];
er = errorbar(y,x,errlow,'.',"horizontal");
er.LineWidth = 1.5;
er.Color = 'Blue';
for i = 1:length(x)
    if p(i) < 0.01
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),"**",'Color','black','FontSize',12)
    elseif  p(i) < 0.05
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),'*','Color','black','FontSize',12)
    end
end
xlim([min(y)-1,max(y)+1])
text(0.6,x(2),"*p < 0.05","FontSize",12)
text(0.6,x(1),"**p < 0.01","FontSize",12)
ctitle = 'Irregular vs Regular Foci Comparison';
title(ctitle)
xlim([-0.2 0.9])
xlabel('Normalized Mean Foci Intensity Ratio')
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.pdf'],'Resolution',600)

%% Bridge Pairing
bridgeinfo = readcell("Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\3D_Bridge_Pair_Coords.xlsx");
% Load all results
% Channel Mean
for i = 1:40
    cmeanfile = [csvfilelist filenamestart '_Intensity_Mean_Ch=' num2str(i) '_Img=1.csv'];
    ctable = readtable(cmeanfile);
    markersname = markers{i};
    markersname(markersname == 'κ') = 'k';
    for ri = 1:length(ctable.IntensityMean)
        Results_bridgeending(ri).(markersname) = ctable.IntensityMean(ri);
    end
end

celltypes = {'AC','MES','OPC','NPC'};
celltypes_markers = {'SOX9','NDRG1','SOX10','ELAVL4'};
celltypeslist = [];
for typei=1:4
    celltypemean(typei) = mean([Results_bridgeending.(celltypes_markers{typei})]);
end
for cellid = 1:length(Results_bridgeending)
    for typei = 1:4
        celltypescore(typei) = Results_bridgeending(cellid).(celltypes_markers{typei}) ...
            ./ celltypemean(typei);
    end
    [~,assigntypei] = max(celltypescore);
    celltypeslist(cellid) = assigntypei; 
end

for axis = {'X','Y','Z'}
    locfile = [csvfilelist filenamestart '_Position_' axis{1} '.csv'];
    csvoutcell = readtable(locfile);
    cpos = table2array([csvoutcell(:,1)]);
    for ri = 1:length([Results_bridgeending.Hoechst1])
        Results_bridgeending(ri).(axis{1}) = cpos(ri);
    end
end

bafposmask = [Results_bridgeending.BANF1]' > 250;
for ri = 1:length([Results_bridgeending.Hoechst1])
    Results_bridgeending(ri).Rupture = bafposmask(ri);
end
xlist = [Results_bridgeending.X];
ylist = [Results_bridgeending.Y];
zlist = [Results_bridgeending.Z];
coordlist = [xlist;ylist;zlist]';
i = 1;
k = [];
list= [];
dist = [];
for bi = 2:length(bridgeinfo)
    ccoord_string(1) = bridgeinfo(bi,7);
    ccoord_string(2) = bridgeinfo(bi,8);
    ccoord = [];
    rem = [];
    if any(~ismissing(ccoord_string{1})) && any(~ismissing(ccoord_string{2}))... 
            && ~contains(ccoord_string(1),'MN') && ~contains(ccoord_string(2),'MN')  
        [ccoord{1,1}, rem{1}] = strtok(ccoord_string{1},'-');
        [ccoord{1,2}, rem{2}] = strtok(ccoord_string{2},'-');
        [ccoord{2,1}, rem{1}] = strtok(rem{1},'-');
        [ccoord{2,2}, rem{2}] = strtok(rem{2},'-');
        [ccoord{3,1}] = strtok(rem{1},'-');
        [ccoord{3,2}] = strtok(rem{2},'-');    
        coordlist_bridgestart = [str2num(ccoord{1,1}) str2num(ccoord{2,1}) str2num(ccoord{3,1})];
        coordlist_bridgeend = [str2num(ccoord{1,2}) str2num(ccoord{2,2}) str2num(ccoord{3,2})];
        [k(i,1),dist(i,1)] = dsearchn(coordlist,coordlist_bridgestart);
        [k(i,2),dist(i,2)] = dsearchn(coordlist,coordlist_bridgeend);
        if dist(i,1) > 10 || dist(i,2) > 10
            k(i,:) = [];
            dist(i,:) = [];
        else
            list{i,1} = celltypes(celltypeslist(k(i,1))); 
            list{i,2} = celltypes(celltypeslist(k(i,2))); 
             i = i +1;
        end
    end
end
% 
% 
% xlist = [];
% ylist = [];
% for i = 1:length(Lineages)
%     clineage = Lineages{i};
%     group = [Results_bridgeending.(clineage)];
%     xlist = [xlist group];
%     ylist = [ylist repmat({clineage}, 1, length(group))];
%     point(i,1) = group(k(1,1))
%     point(i,2) = group(k(1,2))
% end
% boxplot(xlist,ylist)

% Interferon Activation
cmarker = Interferons{1};
list = [];
for i = 1:length(k)
    meanintcell1 = Results_bridgeending(k(i,1)).(cmarker);
    meanintcell2 = Results_bridgeending(k(i,2)).(cmarker);
    list = [list ;[meanintcell1 meanintcell2]];
end

mainmarkers = Interferons; 
errlow = [];
diffperc = [];
for i = 1:length(mainmarkers)
    group = log10([Results_bridgeending.(mainmarkers{i})]./([Results_bridgeending.Hoechst1]).*mean([Results_bridgeending.Hoechst1])+1);% 
    mask = k(:);
    grouppos = group(mask ); %& sox2mask& sox2mask
    groupneg = group;
    [~,p(i),est(i,1:2),~] = ttest2(grouppos,groupneg);
    diffperc(i) = mean(grouppos)/mean(groupneg)-1;
    error = est(i,1:2)/mean(groupneg)+1;
    errlow(i) = mean(grouppos)/mean(groupneg) - error(1);
end
cf = figure;
hold on
y = diffperc;
index = [];
[y,index] = sort(y);
xsorted = mainmarkers;
xsorted = xsorted(index);
x = categorical(xsorted);
x = reordercats(x,xsorted);
errlow = errlow(index);
p = p(index);
barh(x,y);
cf.Position = [1000 200 700 300];
er = errorbar(y,x,errlow,'.',"horizontal");
er.LineWidth = 1.5;
er.Color = 'Blue';
for i = 1:length(x)
    if p(i) < 0.01
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),"**",'Color','black','FontSize',12)
    elseif  p(i) < 0.05
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),'*','Color','black','FontSize',12)
    end
end
xlim([min(y)-1,max(y)+1])
text(0.6,x(2),"*p < 0.05","FontSize",12)
text(0.6,x(1),"**p < 0.01","FontSize",12)
ctitle = 'Difference in Primary Nucleus with Bridge vs without Bridge';
title(ctitle)
xlim([-0.8 0.8])
xlabel('Normalized Mean Nuclear Intensity Ratio')

exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Interferons_on_PN_with_Bridge.pdf'],'Resolution',600)

%% Bridge Pairing All
bridgeinfo = readcell("Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\3D_Bridge_Pair_Coords.xlsx");
% Load all results
% Channel Mean
temp = Results_bridgeending;
for i = 1:40
    cmeanfile = [csvfilelist filenamestart '_Intensity_Mean_Ch=' num2str(i) '_Img=1.csv'];
    ctable = readtable(cmeanfile);
    markersname = markers{i};
    markersname(markersname == 'κ') = 'k';
    for ri = 1:length(ctable.IntensityMean)
        Results_bridgeending(ri).(markersname) = ctable.IntensityMean(ri);
    end
end

celltypes = {'AC','MES','OPC','NPC'};
celltypes_markers = {'SOX9','NDRG1','SOX10','ELAVL4'};
celltypeslist = [];
for typei=1:4
    celltypemean(typei) = mean([Results_bridgeending.(celltypes_markers{typei})]);
end
for cellid = 1:length(Results_bridgeending)
    for typei = 1:4
        celltypescore(typei) = Results_bridgeending(cellid).(celltypes_markers{typei}) ...
            ./ celltypemean(typei);
    end
    [~,assigntypei] = max(celltypescore);
    celltypeslist(cellid) = assigntypei; 
end

for axis = {'X','Y','Z'}
    locfile = [csvfilelist filenamestart '_Position_' axis{1} '.csv'];
    csvoutcell = readtable(locfile);
    cpos = table2array([csvoutcell(:,1)]);
    for ri = 1:length([Results_bridgeending.Hoechst1])
        Results_bridgeending(ri).(axis{1}) = cpos(ri);
    end
end

bafposmask = [Results_bridgeending.BANF1]' > 250;
for ri = 1:length([Results_bridgeending.Hoechst1])
    Results_bridgeending(ri).Rupture = bafposmask(ri);
end
xlist = [Results_bridgeending.X];
ylist = [Results_bridgeending.Y];
zlist = [Results_bridgeending.Z];
coordlist = [xlist;ylist;zlist]';
i = 1;
k = [];
list= [];
dist = [];
for bi = 2:length(bridgeinfo)
    for endi = 1:2
        ccoord_string = bridgeinfo(bi,6+endi);
        ccoord = [];
        rem = [];
        if any(~ismissing(ccoord_string{1})) ... 
                && ~contains(ccoord_string{1},'MN')   
            [ccoord{1,1}, rem{1}] = strtok(ccoord_string{1},'-');  
            [ccoord{2,1}, rem{1}] = strtok(rem{1},'-');     
            [ccoord{3,1}] = strtok(rem{1},'-');          
            coordlist_bridgestart = [str2num(ccoord{1,1}) str2num(ccoord{2,1}) str2num(ccoord{3,1})];
            [k(i,1),dist(i,1)] = dsearchn(coordlist,coordlist_bridgestart);
            if dist(i,1) > 10
                k(i,:) = [];
                dist(i,:) = [];
            else
                list{i,1} = celltypes(celltypeslist(k(i,1))); 
                 i = i +1;
            end
        end
    end
end
% 
% 
% xlist = [];
% ylist = [];
% for i = 1:length(Lineages)
%     clineage = Lineages{i};
%     group = [Results_bridgeending.(clineage)];
%     xlist = [xlist group];
%     ylist = [ylist repmat({clineage}, 1, length(group))];
%     point(i,1) = group(k(1,1))
%     point(i,2) = group(k(1,2))
% end
% boxplot(xlist,ylist)

% Interferon Activation
cmarker = Interferons{1};
list = [];
for i = 1:length(k)
    meanintcell1 = Results_bridgeending(k(i)).(cmarker);
    list = [list ;[meanintcell1]];
end

mainmarkers = [Interferons {'pH2AX'}]; 

baf2mask = [Results_bridgeending.BANF1] > 4000;
errlow = [];
diffperc = [];
for i = 1:length(mainmarkers)
    group = log10([Results_bridgeending.(mainmarkers{i})]./([Results_bridgeending.Hoechst1]).*mean([Results_bridgeending.Hoechst1])+1);% 
    mask = zeros(length(group),1);
    mask(k(:)) = 1;
    mask = mask>0;
    grouppos = group(mask); %& sox2mask& sox2mask
    groupneg = group(~mask & ~baf2mask');
    [~,p(i),est(i,1:2),~] = ttest2(grouppos,groupneg);
    diffperc(i) = mean(grouppos)/mean(groupneg)-1;
    error = est(i,1:2)/mean(groupneg)+1;
    errlow(i) = mean(grouppos)/mean(groupneg) - error(1);
end
cf = figure;
hold on
y = diffperc;
index = [];
[y,index] = sort(y);
xsorted = mainmarkers;
xsorted = xsorted(index);
x = categorical(xsorted);
x = reordercats(x,xsorted);
errlow = errlow(index);
p = p(index);
barh(x,y);
cf.Position = [1000 200 700 300];
er = errorbar(y,x,errlow,'.',"horizontal");
er.LineWidth = 1.5;
er.Color = 'Blue';
for i = 1:length(x)
    if p(i) < 0.01
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),"**",'Color','black','FontSize',12)
    elseif  p(i) < 0.05
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),'*','Color','black','FontSize',12)
    end
end
xlim([min(y)-1,max(y)+1])
text(0.6,x(2),"*p < 0.05","FontSize",12)
text(0.6,x(1),"**p < 0.01","FontSize",12)
ctitle = 'Difference in Primary Nucleus with Bridge vs without Bridge';
title(ctitle)
xlim([-0.8 0.8])
xlabel('Normalized Mean Nuclear Intensity Ratio')

exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Interferons_on_PN_with_Bridge.pdf'],'Resolution',600)



%% Hoechst STD
for i = 1:10
    CVlist(i,:) = [Results_morphology.(['Hoechst' num2str(i) '_std'])]./[Results_morphology.(['Hoechst' num2str(i)])];
end
for i = 1:10
    Hoechstlist(i,:) = [Results_morphology.(['Hoechst' num2str(i) ])];
end
for i = 1:10
    STDlist(i,:) = [Results_morphology.(['Hoechst' num2str(i) '_std'])];
end

figure,
testcellrandom = [343,143,323,221,144,526,451,145,233,31];
testchannel = 1:10;
hold on
for i = testcellrandom
    x = [testchannel];
    y = CVlist(testchannel,i);
    plot(x,y)
end
ylabel('CV')
xlabel('Cycle')
ctitle = 'CV of 10 nucleus over 10 cycles';
title(ctitle)
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.jpg'],'Resolution',600)
hold off


figure,
hold on
for i = testcellrandom
    x = [testchannel];
    y = Hoechstlist(testchannel,i);
    plot(x,y)
end
ylabel('Mean Hoechst Intensity')
xlabel('Cycle')
ctitle = 'Mean Hoechst of 10 nucleus over 10 cycles';
title(ctitle)
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.jpg'],'Resolution',600)
hold off

figure,
hold on
for i = testcellrandom
    x = [testchannel];
    y = STDlist(testchannel,i);
    plot(x,y)
end
ylabel('Std Hoechst Intensity')
xlabel('Cycle')
ctitle = 'Std of 10 nucleus over 10 cycles';
title(ctitle)
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' ctitle '.jpg'],'Resolution',600)
hold off

%% DNA comparison
Metrics = {'Volume','Sphericity','Ellipticityoblate','Ellipticityprolate','IrregularIndex'};
celltypes = {'AC','MES','OPC','NPC'};
celltypes_markers = {'SOX9','NDRG1','SOX10','ELAVL4'};
celltypes_markers_cutoff = [371,600,800,200];
sox2mask = [Results_morphology.SOX2mask]' == 1;
for typei=1:4
    celltypemean(typei) = mean([Results_morphology.(celltypes_markers{typei})]);
end
for cellid = 1:length(Results_morphology)
    for typei = 1:4
        celltypescore(typei) = Results_morphology(cellid).(celltypes_markers{typei}) ...
            ./ celltypes_markers_cutoff(typei);
%         if typei == 2
%             mean1 = Results_morphology(cellid).('NDRG1') ./ mean([Results_morphology.('NDRG1')]);
%             mean2 = Results_morphology(cellid).('HIF1a') ./ mean([Results_morphology.('HIF1a')]);
%             celltypescore(typei) = mean([mean1,mean2]);
%         end
    end
    [~,assigntypei] = max(celltypescore);
    celltypeslist(cellid) = assigntypei; 
%     for typei = 1:4
%         if Results_morphology(cellid).(celltypes_markers{typei}) > celltypes_markers_cutoff(typei)
%             celltypeslist(cellid) = typei;
%         end
%     end
    [~,~,ix] = unique(celltypeslist);
    C = accumarray(ix,1).';
end

for i = 1:3
    
    f = figure;
    f.Position = [1000 500 400 400];
    switch i
        case 1
            cmetric = 'Hoecst1 Mean';
        case 2
            cmetric = 'Hoecst1 Standard Deviation';
        case 3
            cmetric = 'Hoecst1 CV';
    end
    hold on
    data = [];
    cat = [];
    for typei = 1:4
        cmask = celltypeslist' == typei;
        switch i
            case 1
                cdata = [Results_morphology.Hoechst1];
            case 2
                cdata = [Results_morphology.Hoechst1_std];
            case 3
                cdata = [Results_morphology.Hoechst1_std]./[Results_morphology.Hoechst1];
        end
        posmask = cdata'>0;
        cdata_pos =  cdata(cmask& sox2mask & posmask);
        data = [data cdata_pos ];
        violindata{typei} = cdata_pos; 
        cat = [cat repmat({[celltypes{typei}]},1,length(cdata_pos))];

        grouppos = cdata_pos;
        groupneg = cdata(~cmask& sox2mask & posmask);
        [~,p,est,~] = ttest2(grouppos,groupneg);
        diffperc = abs(mean(grouppos)/mean(groupneg)*100-100);
        if p < 0.1
            if mean(est) > 0
                text(typei+0.1,mean(grouppos)*1.2,['p = ' num2str(p) newline ... 
                    num2str(diffperc,2) '% Higher' newline ...
                    'in ' celltypes{typei}],'Color','black','FontSize',9)
            else
                text(typei+0.1,mean(grouppos)*1.2,['p = ' num2str(p) newline ...
                    num2str(diffperc,2) '% Lower' newline ...
                    'in ' celltypes{typei}],'Color','black','FontSize',9)
            end
        end
    end
    %boxplot(data,cat)
    xlabel = celltypes;
    violin(violindata,'xlabel',xlabel)
    
    title(cmetric)
    hold off
    
    exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/' cmetric '.jpg'],'Resolution',600)
end

%% Irregularity Index vs Lamin
f = figure;
f.Position = [500 500 400 200];
xlim([2 4])
clear xlabel ylabel
metric = 'Irregularity Index';
marker = 'LaminB1';
x = [Results_morphology.Volume]./[Results_morphology.Sphericity];
%y = [Results_morphology.(marker)]./[Results_morphology.Area].*mean([Results_morphology.Area]);
y = [Results_morphology.(marker)];%.*[Results_morphology.Volume]./[Results_morphology.Area];
sox2mask = [Results_morphology.SOX2mask];
x = x(y>50&sox2mask);
y = y(y>50&sox2mask);
y = log10(y);
% scatter(x,y)
% corrcoef(x,y)
% xlabel('Irregular Index')
% ylabel(marker)

mask_low = x < prctile(x,50);
mask_high = x > prctile(x,50);
y_low = y(mask_low);
y_high = y(mask_high);
mean(y_low)
mean(y_high)
[~,p]= ttest2(y_low,y_high)
hold on
ksdensity(y_low)
ksdensity(y_high)
xlabel(marker)
ylabel('Probabilities')
legend({'Low AI','High AI'})
hold off
%boxplot([x_low' x_high'])
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images\AIvs_' marker '.jpg'],'Resolution',600)

f = figure;
f.Position = [500 500 400 200];
xlim([0 5])
marker = 'LAMINAC';
x = [Results_morphology.Volume]./[Results_morphology.Sphericity];
%y = [Results_morphology.(marker)];%./[Results_morphology.Area].*mean([Results_morphology.Area]);
y = [Results_morphology.(marker)];%.*[Results_morphology.Volume]./[Results_morphology.Area];
sox2mask = [Results_morphology.SOX2mask]>=0;
x = x(y>50&sox2mask);
y = y(y>50&sox2mask);
y = log10(y);
% scatter(x,y)
% corrcoef(x,y)
% xlabel('Irregular Index')
% ylabel(marker)

mask_low = x < prctile(x,20);
mask_high = x > prctile(x,80);
y_low = y(mask_low);
y_high = y(mask_high);
mean(y_low)
mean(y_high)
[~,p]= ttest2(y_low,y_high)
hold on
ksdensity(y_low)
ksdensity(y_high)
xlabel(marker)
ylabel('Probabilities')
legend({'Low AI','High AI'})
hold off
%boxplot([x_low' x_high'])
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images\AIvs_' marker '.jpg'],'Resolution',600)

%% Sphericity vs Volume
f = figure;
f.Position=[500 500 200 200];
hold on
xmetric = 'Sphericity';
ymetric = 'Volume';
rupturemask = [Results_morphology.Rupture];
x = [Results_morphology.(xmetric)];
y = [Results_morphology.(ymetric)];
xrup = x(rupturemask);
yrup = y(rupturemask);
scatter(xrup,yrup,7,'red','filled')
%densityplot(x,y,[20 20])
xnrup = x(~rupturemask);
ynrup = y(~rupturemask);
scatter(xnrup,ynrup,7,'blue','filled')
%densityplot2(x,y,[20 20])
%alpha(.5)
xlabel(xmetric)
ylabel(ymetric)
corrcoef(x,y)
hold off
 exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/SphericityvsVolume_PN.pdf'],'Resolution',600)

f = figure;
f.Position=[500 500 200 200];
hold on
xmetric = 'Sphericity';
ymetric = 'Volume';
rupturemask = [Results_micro.BANF1] > 300;
x = [Results_micro.(xmetric)];
y = [Results_micro.(ymetric)];
xrup = x(rupturemask);
yrup = y(rupturemask);
scatter(xrup,yrup,4,'red','filled')
xnrup = x(~rupturemask);
ynrup = y(~rupturemask);
scatter(xnrup,ynrup,4,'blue','filled')
xlabel(xmetric)
ylabel(ymetric)
corrcoef(x,y)
hold off
 exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/SphericityvsVolume_MN.pdf'],'Resolution',600)

f = figure;
f.Position=[500 500 400 400];
xmetric = 'Sphericity';
ymetric = 'Volume';
rupturemask = [Results_morphology.Rupture];
x = [Results_morphology.(xmetric)];
y = [Results_morphology.(ymetric)];
xrup = x(rupturemask);
yrup = y(rupturemask);
xgrid=linspace(0.65,1);
ygrid=linspace(0,800);
[x1,y1] = meshgrid(xgrid, ygrid);
% Perform kernel density estimate
% [x y] is actual data, xi is the desired grid points to evaluate
% f is an estimate of the density, ep(:,1) is the X location, ep(:,2) is the y location
xi = [x1(:) y1(:)];
[f,ep]=ksdensity([xrup' yrup'],xi); % remove the outputs to see a 3D plot of the distribution
% format data in matrix for contourf and plot
X = reshape(ep(:,1),length(xgrid),length(ygrid));
Y = reshape(ep(:,2),length(xgrid),length(ygrid));
Z = reshape(f,length(xgrid),length(ygrid));
contour(X,Y,Z,10)
colormap(slanCM('Blues'))
clim([-0.022 0.022])
ax=gca;
ax.XAxis.Exponent = 0;
xtickformat('%.4f')
xlabel(xmetric)
ylabel(ymetric)
colorbar
 exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/SphericityvsVolume_PN_contour_BAF+.pdf'],'Resolution',600)

f = figure;
f.Position=[500 500 400 400];
xnrup = x(~rupturemask);
ynrup = y(~rupturemask);
 % Define grid
xgrid=linspace(0.65,1);
ygrid=linspace(0,800);
[x1,y1] = meshgrid(xgrid, ygrid);
% Perform kernel density estimate
% [x y] is actual data, xi is the desired grid points to evaluate
% f is an estimate of the density, ep(:,1) is the X location, ep(:,2) is the y location
xi = [x1(:) y1(:)];
[f,ep]=ksdensity([xnrup' ynrup'],xi); % remove the outputs to see a 3D plot of the distribution
% format data in matrix for contourf and plot
X = reshape(ep(:,1),length(xgrid),length(ygrid));
Y = reshape(ep(:,2),length(xgrid),length(ygrid));
Z = reshape(f,length(xgrid),length(ygrid));
contour(X,Y,Z,10)
colormap(slanCM('Greys'))
clim([-0.022 0.022])
ax=gca;
ax.XAxis.Exponent = 0;
xtickformat('%.4f')
xlabel(xmetric)
ylabel(ymetric)
corrcoef(x,y)
colorbar
 exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/SphericityvsVolume_PN_contour_BAF-.pdf'],'Resolution',600)

f = figure;
f.Position=[500 500 400 400];
xmetric = 'Sphericity';
ymetric = 'Volume';
rupturemask = [Results_micro.BANF1] > 300;
x = [Results_micro.(xmetric)];
y = [Results_micro.(ymetric)];
xrup = x(rupturemask);
yrup = y(rupturemask);
xgrid=linspace(0.65,1);
ygrid=linspace(0,25);
[x1,y1] = meshgrid(xgrid, ygrid);
% Perform kernel density estimate
% [x y] is actual data, xi is the desired grid points to evaluate
% f is an estimate of the density, ep(:,1) is the X location, ep(:,2) is the y location
xi = [x1(:) y1(:)];
[f,ep]=ksdensity([xrup' yrup'],xi); % remove the outputs to see a 3D plot of the distribution
% format data in matrix for contourf and plot
X = reshape(ep(:,1),length(xgrid),length(ygrid));
Y = reshape(ep(:,2),length(xgrid),length(ygrid));
Z = reshape(f,length(xgrid),length(ygrid));
contour(X,Y,Z,10)
colormap(slanCM('Oranges'))
clim([-0.5 0.9])
ax=gca;
ax.XAxis.Exponent = 0;
xtickformat('%.4f')
xlabel(xmetric)
ylabel(ymetric)
colorbar
 exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/SphericityvsVolume_MN_contour_BAF+.pdf'],'Resolution',600)

f = figure;
f.Position=[500 500 400 400];
xnrup = x(~rupturemask);
ynrup = y(~rupturemask);
 % Define grid
xgrid=linspace(0.65,1);
ygrid=linspace(0,20);
[x1,y1] = meshgrid(xgrid, ygrid);
% Perform kernel density estimate
% [x y] is actual data, xi is the desired grid points to evaluate
% f is an estimate of the density, ep(:,1) is the X location, ep(:,2) is the y location
xi = [x1(:) y1(:)];
[f,ep]=ksdensity([xnrup' ynrup'],xi); % remove the outputs to see a 3D plot of the distribution
% format data in matrix for contourf and plot
X = reshape(ep(:,1),length(xgrid),length(ygrid));
Y = reshape(ep(:,2),length(xgrid),length(ygrid));
Z = reshape(f,length(xgrid),length(ygrid));
contour(X,Y,Z,10)
colormap(slanCM('Greys'))
clim([-6 6])
ax=gca;
ax.XAxis.Exponent = 0;
xtickformat('%.4f')
xlabel(xmetric)
ylabel(ymetric)
colorbar
 exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/SphericityvsVolume_MN_contour_BAF-.pdf'],'Resolution',600)

 %% 3D Lamin at Rupture site

 % Map Foci to Primary nucleus
xlist = [Results_all.X];
ylist = [Results_all.Y];
zlist = [Results_all.Z];
coordlist = [xlist;ylist;zlist]';
xlist_Foci = [Results_Foci.X];
ylist_Foci = [Results_Foci.Y];
zlist_Foci = [Results_Foci.Z];
coordlist_Foci = [xlist_Foci;ylist_Foci;zlist_Foci]';
[k_Foci,dist] = dsearchn(coordlist,coordlist_Foci);

marker = 'LAMINAC';
Group_foci = [Results_Foci.(marker)];%./[Results_Foci.Hoechst1];
Group_Nuc = [Results_all.(marker)];%./[Results_all.Hoechst1];
Group_Nuc = Group_Nuc(k_Foci);
ratio = Group_foci./Group_Nuc;
ratio(~(ratio>0)) = [];
mean(ratio)
%% Temp
x = [Results_morphology.Volume];
x = x(x>2);
ksdensity(x,'support','positive')
exportgraphics(gcf,['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_05_GBM_3D\Nuclear_Shape\Images/Volume_PN_density.pdf'],'Resolution',600)
