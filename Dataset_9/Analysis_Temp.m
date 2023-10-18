% Processing Multiple slides
clear all
filename.tissues = {'LSP15550', 'LSP15558'...
                    };    

for ctissuename_index = 1:length(filename.tissues)
    % 1) THE FOLDER DIRECTORIES
    ctissuename = filename.tissues{ctissuename_index}; % Current Tissue Name
    filename.tissue = ctissuename; 
    filename.datafolder = ['Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Liposarcoma\' ctissuename '\']; %Main Analysis location
    filename.rawfolder =  ['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Liposarcoma\' ctissuename '\']; %Ometif raw file location
    filename.graphsfolder = 'Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Liposarcoma\Graphs\'; % Output Graph Folder
    filename.suffix = '.tif';
    filename.primaryNucfolder = [filename.datafolder 'ANALYSIS20230830\']; % Primary Nucleus Analysis folder
    filename.microNucfolder = [filename.datafolder 'ANALYSIS20230830micro\']; % Micronucleus Analysis folder
    filename.resultfolder = [filename.datafolder 'Results\']; % Result folder
    filename.fullstackfolder = [filename.datafolder 'FullStacks' filesep]; % Full Stacks Folder
    filename.folders.fullstacks = 'FullStacks\';
    filename.folders.coordinates = 'Coordinates_FullStacks\';
    filename.dim = 2; % Data Dimension
    
    % 2) SLIDE SPECIFIC PARAMETERS:
    if exist(filename.fullstackfolder,"dir")
        tempdir = dir(filename.fullstackfolder);
        [~,ind]=sort([tempdir.datenum]);
        tempdir = tempdir(ind);
        coren = 0;
        for i1 = 1:(length(dir(filename.fullstackfolder)))
            if contains(tempdir(i1).name,filename.tissues{ctissuename_index})
                coren = coren + 1;
                filename.realcoreinfo(coren).name = tempdir(i1).name;
                filename.realcoreinfo(coren).index = coren;
            end
        end
        clear i1 tempdir coren ind
    end
    filename.cycles = 12;  % cycles to analyse
    filename.maxround = filename.cycles*4;
    filename.ilastiksegfol = 'IlastikSegmentation\'; % Ilastik Folder setup
    filename.ilastiksuffix = '_Probabilities.tif';
    filename.segsuffix = '_Seg.tif';                    
    filename.ilastikround = 12; % Rounds to be used to segment in Ilastik (<= maxround)
    
    % 3) USER DESIRED PARAMETERS 
    filename.sizefield = 5000; %Size of field desired for a square 
    filename.crops = 1; %# of cropped fields desired per field 
    options.crops.DAPIbkgd_crops = 100;
    options.crops.prctile = 75;
    options.DAPIbkgd = 100; % thresh for DAPI used to check if we keep the field
    options.DAPIbkgd_crops75prctile = options.DAPIbkgd*5;
    
    % 4) OPTIONS  
    % Step 3: SEGMENTATION OPTIONS 
    options.nuc = 1;
    options.cyt = 2;
    options.backgr = 3;
    options.max_prob = 65535;
    options.bkgdprob_min = 0.2;
    options.cytprob_min  = 0.9;
    options.nucprob_min  = 0.7;
    options.nucprob_min_micro  = 0.6;
    options.pausetime = 0;
    
    
    % Step 4: MEASUREMENTS AND FOCI OPTIONS 
    options.tilesize = 1024;
    options.scale = 8;
    options.buffer = 256;
    options.coreradius = 200;
    options.dim = ['%02d'];
    options.rows = 2;
    options.cols = 2;
    options.background = 500;
    options.sigma = 0.3;
    options.alpha = 0.1;
    options.fileordchange = [1 2 3 4];
    options.cellsize = 100;
    options.cytosize = options.cellsize/5;
    options.corrthresh = 0.5;
    options.appendflag = 0;
    options.dates = '20230830_';
    
    filename.resultfile = [options.dates 'Results.mat'];
    filename.resultfilemicro = [options.dates 'microResults.mat'];
    
    
    options.substructures.flag = 1;
    options.substructures.channel = 20;
    options.substructures.bkgsize = 15;
    options.substructures.SV_thresh = 1;
    options.substructures.nucdilation = 3; % in case the nuclear segmentation is very stringent
    
    % 5) Cylinter
    filename.cylinterfol = 'Cylinter\';
    filename.cylinterinputfol = 'Input\';
    filename.cylintermaskfol = 'mask';
    filename.cylintersegfol = 'seg';
    filename.cylintertiffol = 'tif';
    
    % 6) Roi from Omero
    filename.roifile = [filename.datafolder 'Batch_ROI_Export'];
    
    % 7) Aggregate DATA
    options.cgaschannel = 42;

    load([filename.resultfolder 'Aggregate_' filename.resultfile])
    temp{ctissuename_index} = AggrResults_primary; 
    temp_MN{ctissuename_index} = AggrResults_micro; 
end
temp{1} = rmfield(temp{1},'Corenum');
Fname =  fieldnames(temp{1});
AggrResults_primary = [];
for f = 1:length(Fname)
    AggrResults_primary.(Fname{f}) = [];
    for i = 1:length(temp)
        if ~isempty(temp{i}) && ~isempty(temp{i}.(Fname{f}))
            [mr,mc] = size(temp{i}.(Fname{f}));
            if (mr > 1)
                AggrResults_primary.(Fname{f}) = [AggrResults_primary.(Fname{f});temp{i}.(Fname{f})];
            end
        end
    end
end
Fname =  fieldnames(temp_MN{1});
AggrResults_micro = [];
for f = 1:length(Fname)
    AggrResults_micro.(Fname{f}) = [];
    for i = 1:length(temp_MN)
        if ~isempty(temp_MN{i}) && ~isempty(temp_MN{i}.(Fname{f}))
            [mr,mc] = size(temp_MN{i}.(Fname{f}));
            if (mr > 1)
                AggrResults_micro.(Fname{f}) = [AggrResults_micro.(Fname{f});temp_MN{i}.(Fname{f})];
            end
        end
    end
end
Channelcutoffcell = readcell([filename.datafolder 'Channelscutoff_BC.xlsx']);
[maxchannel,~] = size(Channelcutoffcell);
AllChannel = Channelcutoffcell(:,1);
AllChanneln = 1:maxchannel;
Channelcutoff = cell2mat(Channelcutoffcell(:,2));
amps = cell2mat(Channelcutoffcell(:,3));
%Channelchar = Channelcutoffcell(:,4);
ChannelL = Channelcutoffcell(:,5);

mdm2mask = AggrResults_primary.MeanNucSign(:,6) > 400;
%mdm2mask = AggrResults_primary.MeanNucSign(:,6) > 400 & AggrResults_primary.MeanNucSign(:,6) > 2* AggrResults_primary.MeanCytSign(:,6);

%% Primary Nucleus with Rutpure versus no Rupture
clear resulttable
ChannelTypes = {'BAFfriends','Interferon','Other'};
for ct = 1:length(ChannelTypes)
    ChannelType = ChannelTypes{ct};
    c = [];
    for i = 1:length(ChannelL)
        if length(ChannelL{i})>1 && strcmp(ChannelL{i},ChannelType)
            c = [c i];
        end
    end
    % c([2]) = [];
    resulttable = [];
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        caggr = AggrResults_primary.MeanFullCellSign(:,ccn);%./AggrResults_primary.MeanNucSign(:,1);
        posmask = caggr>0;
        cposmask = AggrResults_primary.AreaSubstruct(:,1) > 0;
        cnegmask = AggrResults_primary.AreaSubstruct(:,1) == 0;
        cmeanBAFpos = caggr(cposmask&posmask&mdm2mask);
        cmeanBAFneg = caggr(cnegmask&posmask&mdm2mask);
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
 
    figure,
    hold on
    y = [resulttable{2,2:end}]-1;
    [y,index] = sort(y);
    xsorted = resulttable(1,2:end);
    xsorted = xsorted(index);
    x = categorical(xsorted);
    x = reordercats(x,xsorted);
    errlow = errlow(index);
    barh(x,y);
    er = errorbar(y,x,errlow,'.',"horizontal");
    er.LineWidth = 1.5;
    er.Color = 'Blue';
    for i = 1:length(x)
        if resulttable{3,index(i)+1} < 0.01
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),"**",'Color','black','FontSize',12)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',12)
        end
    end
    xlim([min(y)-0.2,max(y)+0.25])
    text(0.1,x(2),"*p < 0.05","FontSize",12)
    text(0.1,x(1),"**p < 0.01","FontSize",12)
    title(['Difference in Primary Nucleus with Rutpure versus no Rupture'])
    xlabel('Normalized Mean Nuclear Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    exportgraphics(gca,[filename.graphsfolder 'Difference in Primary Nucleus with Rutpure versus no Rupture in ', ChannelType,'.jpg'],'Resolution',600,'ContentType','vector')
end

%% Primary Nucleus with Micronucleus versus with no Micronucleus
clear resulttable
ChannelTypes = {'BAFfriends','Interferon','Other'};
for ct = 1:length(ChannelTypes)
    ChannelType = ChannelTypes{ct};
    c = [];
    for i = 1:length(ChannelL)
        if length(ChannelL{i})>1 && strcmp(ChannelL{i},ChannelType)
            c = [c i];
        end
    end
    resulttable = [];
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        caggr = AggrResults_primary.MeanNucSign(:,ccn);%./AggrResults_primary.MeanNucSign(:,1);
        posmask = caggr>0;
        cposmask = AggrResults_primary.micronucleus(:,1)>0;
        cnegmask = AggrResults_primary.micronucleus(:,1)==0;
        cmeanBAFpos = caggr(cposmask&posmask&mdm2mask);
        cmeanBAFneg = caggr(cnegmask&posmask&mdm2mask);
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    figure,
    hold on
    y = [resulttable{2,2:end}]-1;
    [y,index] = sort(y);
    xsorted = resulttable(1,2:end);
    xsorted = xsorted(index);
    x = categorical(xsorted);
    x = reordercats(x,xsorted);
    errlow = errlow(index);
    barh(x,y);
    er = errorbar(y,x,errlow,'.',"horizontal");
    er.LineWidth = 1.5;
    er.Color = 'Blue';
    for i = 1:length(x)
        if resulttable{3,index(i)+1} < 0.01
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),"**",'Color','black','FontSize',12)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',12)
        end
    end
    xlim([min(y)-0.2,max(y)+0.25])
    text(0.1,x(2),"*p < 0.05","FontSize",12)
    text(0.1,x(1),"**p < 0.01","FontSize",12)
    title(['Difference in Primary Nucleus with Micronucleus versus with no Micronucleus'])
    xlabel('Normalized Mean Nuclear Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    exportgraphics(gca,[filename.graphsfolder 'Difference in Primary Nucleus with Micronucleus versus with no Micronucleus in ', ChannelType,'.jpg'],'Resolution',600)
end
%% Primary Nucleus with Micronucleus rupture versus no Micronucleus rupture
clear resulttable
ChannelTypes = {'BAFfriends','Interferon','Other'};
for ct = 1:length(ChannelTypes)
    ChannelType = ChannelTypes{ct};
    c = [];
    for i = 1:length(ChannelL)
        if length(ChannelL{i})>1 && strcmp(ChannelL{i},ChannelType)
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    %sox2mask = AggrResults_primary.MeanNucSign(:,15) > 80;
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        aggrcurrent = AggrResults_primary.MeanFullCellSign(:,ccn);
        cposmask = AggrResults_primary.micronucleusrupture==1;
        cnegmask = ~cposmask & AggrResults_primary.micronucleus==1;%& cposmask_cores_mask;
        cmeanBAFpos = aggrcurrent(cposmask&mdm2mask);
        cmeanBAFneg = aggrcurrent(cnegmask&mdm2mask);
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    figure,
    hold on
    y = [resulttable{2,2:end}]-1;
    [y,index] = sort(y);
    xsorted = resulttable(1,2:end);
    xsorted = xsorted(index);
    x = categorical(xsorted);
    x = reordercats(x,xsorted);
    errlow = errlow(index);
    barh(x,y);
    er = errorbar(y,x,errlow,'.',"horizontal");
    er.LineWidth = 1.5;
    er.Color = 'Blue';
    for i = 1:length(x)
        if resulttable{3,index(i)+1} < 0.01
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),"**",'Color','black','FontSize',12)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',12)
        end
    end
    xlim([min(y)-0.2,max(y)+0.25])
    text(0.1,x(2),"*p < 0.05","FontSize",12)
    text(0.1,x(1),"**p < 0.01","FontSize",12)
    title('Difference in Primary Nucleus with Micronucleus rupture versus no Micronucleus rupture')
    xlabel('Normalized Mean Nuclear and Cytoplasmic Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    %exportgraphics(gca,'test.jpg','Resolution',600)
    exportgraphics(gcf,[filename.graphsfolder 'Difference in Primary Nucleus with Micronucleus rupture versus no Micronucleus rupture ', ChannelType,'.jpg'],'Resolution',600,'ContentType','vector');
end

%% Primary Nucleus with Rupture versus no Rupture
clear resulttable
c = [];
for i = 1:length(ChannelL)
    if ~ismissing(ChannelL{i}) & Channelcutoff(i) > 0 & (strcmp(ChannelL(i),'Lineage') ) %InterferonLineage
        c = [c i];
    end
end
resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
resulttable{3,1} = 't2test-pvalue';
%sox2mask = AggrResults_primary.MeanNucSign(:,15) > 80;
% = AggrResults_primary.tumor == 1;
for cn = 1:length(c)
    channeln = c(cn);
    ccn = AllChanneln(channeln);
    resulttable{1,cn+1} = AllChannel{channeln};
%     if Channelchar{channeln} == 'C'
%       aggrcurrent = AggrResults_primary.MeanCytSign(:,ccn);
%     else
%       aggrcurrent = AggrResults_primary.MeanNucSign(:,ccn);
%     end
    aggrcurrent = AggrResults_primary.MeanFullCellSign(:,ccn);
    cposmask = AggrResults_primary.AreaSubstruct(:,1) > 0;
    cnegmask = ~cposmask ;%& cposmask_cores_mask;
    posmask = aggrcurrent > 0;
    cmeanBAFpos = aggrcurrent(cposmask&sox2mask&posmask);
    cmeanBAFneg = aggrcurrent(cnegmask&sox2mask&posmask);
    resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
    [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
    error = ci/mean(cmeanBAFneg)+1;
    errlow(cn) = resulttable{2,cn+1} - error(1);
    errhigh(cn) = error(2) - resulttable{2,cn+1};
    resulttable{3,cn+1} = p;
end

figure,
hold on
y = [resulttable{2,2:end}]-1;
[y,index] = sort(y);
xsorted = resulttable(1,2:end);
xsorted = xsorted(index);
x = categorical(xsorted);
x = reordercats(x,xsorted);
errlow = errlow(index);
barh(x,y);
er = errorbar(y,x,errlow,'.',"horizontal");
er.LineWidth = 1.5;
er.Color = 'Blue';
for i = 1:length(x)
    if resulttable{3,index(i)+1} < 0.01
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),"**",'Color','black','FontSize',12)
    elseif  resulttable{3,index(i)+1} < 0.05
        if y(i) > 0
            plotpoint = y(i)+errlow(i)+0.01;
        else 
            plotpoint = y(i)-errlow(i)-0.04;
        end
        text(plotpoint,x(i),'*','Color','black','FontSize',12)
    end
end
xlim([min(y)-0.2,max(y)+0.25])
text(0.1,x(2),"*p < 0.05","FontSize",12)
text(0.1,x(1),"**p < 0.01","FontSize",12)
title('Difference in Primary Nucleus with Rupture versus no Rupture')
xlabel('Normalized Mean Nuclear and Cytoplasmic Intensity Ratio')
hold off
set(gcf,'position',[400,-10,1000,600])
%exportgraphics(gca,'test.jpg','Resolution',600)
exportgraphics(gcf,'Difference in Primary Nucleus with Rupture versus no Rupture.pdf','Resolution',600,'ContentType','vector')

%%
Bafposcount = sum(AggrResults_primary.AreaSubstruct(:,1)>0);
Bafnegcount = sum(AggrResults_primary.AreaSubstruct(:,1)==0);
%cgas
markposcount(1,1) = sum(AggrResults_primary.MeanFociSign(:,22) > 200 &  (AggrResults_primary.MeanFociSign(:,22) > 1.3*AggrResults_primary.MeanNucSign(:,22)) & AggrResults_primary.AreaSubstruct(:,1)>0);
markposcount(2,1) = sum(AggrResults_primary.MeanFociSign(:,14) > 200 &  (AggrResults_primary.MeanFociSign(:,14) > 1.3*AggrResults_primary.MeanNucSign(:,14)) & AggrResults_primary.AreaSubstruct(:,1)>0);
markposcount(3,1) = sum(AggrResults_primary.MeanNucSign(:,35) > 200 &  (AggrResults_primary.MeanNucSign(:,35) > 1.3*AggrResults_primary.MeanCytSign(:,35)) & AggrResults_primary.AreaSubstruct(:,1)>0);
markposcount(3,2) = sum(AggrResults_primary.MeanNucSign(:,35) > 200 &  (AggrResults_primary.MeanNucSign(:,35) > 1.3*AggrResults_primary.MeanCytSign(:,35)) & AggrResults_primary.AreaSubstruct(:,1)==0);

bafmask = AggrResults_micro.MeanNucSign(:,20) > 200 &  (AggrResults_micro.MeanNucSign(:,20) > 1.3*AggrResults_micro.MeanCytSign(:,20));
Bafposcount_MN = sum(bafmask);
Bafnegcount_MN = sum(~bafmask);
%cgas
markposcount_MN(1,1) = sum(AggrResults_micro.MeanNucSign(:,22) > 200 &  (AggrResults_micro.MeanNucSign(:,22) > 1.3*AggrResults_micro.MeanCytSign(:,22)) & bafmask);
markposcount_MN(2,1) = sum(AggrResults_micro.MeanNucSign(:,14) > 200 &  (AggrResults_micro.MeanNucSign(:,14) > 1.3*AggrResults_micro.MeanCytSign(:,14)) & bafmask);
markposcount_MN(1,2) = sum(AggrResults_micro.MeanNucSign(:,22) > 200 &  (AggrResults_micro.MeanNucSign(:,22) > 1.3*AggrResults_micro.MeanCytSign(:,22)) & ~bafmask);
markposcount_MN(2,2) = sum(AggrResults_micro.MeanNucSign(:,14) > 200 &  (AggrResults_micro.MeanNucSign(:,14) > 1.3*AggrResults_micro.MeanCytSign(:,14)) & ~bafmask);

%% Cross Correlation

c = [];
for i = 1:length(ChannelL)
    if ~ismissing(ChannelL{i}) 
        c = [c i];
    end
end
coef = [];
soxmask = AggrResults_primary.MeanFullCellSign(:,15) > 50;
for i = 1:length(c)
    for j = 1:length(c)
        Ai = c(i);
        Bi = c(j);
        A = AggrResults_primary.MeanNucSign(:,Ai);
        B = AggrResults_primary.MeanNucSign(:,Bi);
        clearnegmask = A>0 & B>0 & soxmask;
        A = A(clearnegmask);
        B = B(clearnegmask);
        A = log10(A)/mean(log10(A));
        B = log10(B)/mean(log10(B));
        matrixab = [A,B];
        ccoef = corrcoef(matrixab);
        coef{i+1,j+1} = round(ccoef(2),2,"significant");
    end
end
coef(2:2+length(c)-1,1) = AllChannel(c);
coef(1,2:2+length(c)-1) = AllChannel(c);

clf
eucD = pdist(cell2mat(coef(2:end,2:end)),'euclidean');
clustTreeEuc = linkage(eucD,'average');
[h,nodes,outperm] = dendrogram(clustTreeEuc,0);

c = c(outperm);
coef = [];
for i = 1:length(c)
    for j = 1:length(c)
        Ai = c(i);
        Bi = c(j);
        A = AggrResults_primary.MeanNucSign(:,Ai);
        B = AggrResults_primary.MeanNucSign(:,Bi);
        clearnegmask = A>0 & B>0 & soxmask ;% &nonimmune;
        A = A(clearnegmask);
        B = B(clearnegmask);
        A = log10(A)/mean(log10(A));
        B = log10(B)/mean(log10(B));
        matrixab = [A,B];
        ccoef = corrcoef(matrixab);
        coef{i+1,j+1} = round(ccoef(2),2,"significant");
    end
end
coef(2:2+length(c)-1,1) = AllChannel(c);
coef(1,2:2+length(c)-1) = AllChannel(c);

%normalize
cmat = cell2mat(coef(2:end,2:end));
%Heatmap
h = heatmap(coef(1,2:end),coef(2:end,1)',cell2mat(coef(2:end,2:end)));
set(gcf,'position',[400,100,1000,800])
colormap('Hot')
exportgraphics(gcf,['Heatmap.jpg'],'Resolution',600)
