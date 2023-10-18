function Dataset_2_Analysis_Bafanalysisaggregate
    clear all
    driveprefix = 'Y:\'; %HiTS drive prefix;
    for foldinginitializationsteps=1
        %Change All parameters below  
        filename.suffix = '.tif';
        filename.datafolder = 'Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Inteferons\';
        filename.analfolder = [filename.datafolder 'ANALYSIS03302023micro\'];
        filename.resultfolder = [filename.analfolder 'Results\'];
        filename.codesfolder = [filename.analfolder 'Codes Used\'];
        filename.fullstackfolder = [filename.datafolder 'dearray' filesep]; 
        
        filename.cycles = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];  % cycles to analyse
        filename.maxcycle = 14;
        
        filename.ilastiksegfol = 'IlastikSegmentation\';
        filename.ilastiksuffix = '_Probabilities.tif';
        filename.segsuffix = '_Seg.tif';
        
        options.nuc = 1;
        options.cyt = 2;
        options.backgr = 3;
        options.max_prob = 65535;
        options.bkgdprob_min = 0.2;
        options.cytprob_min  = 0.8;
        options.nucprob_min  = 0.8;
        options.pausetime = 0;
        
        for i1 = 1:(length(dir(filename.fullstackfolder))-3)
            filename.realcoreinfo(i1).index = i1;
        end
        
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
        options.cellsize = 50;
        options.cytosize = options.cellsize/5;
        options.corrthresh = 0.5;
        options.appendflag = 0;
        options.dates = '2023-4-10_';
        
        filename.resultfile = [options.dates 'Results_PostCylinter.mat'];
        
        options.substructures.flag = 1;
        options.substructures.channel = 14;
        options.substructures.bkgsize = 15;
        options.substructures.SV_thresh = 0.8;
        options.substructures.nucdilation = 3; % in case the nuclear segmentation is very stringent
        
        DAPIslice = (1:filename.maxcycle)*4 -3;
        
            
        %For 10X Montage
        % rows = 1;
        % cols = 1;
        % maxtile = rows*cols;
        cut = 100;
        background = 1000;
        radius10X = 700;
        
        % images per each core
        rows2 = 2;
        cols2 = 2;
        maxtile2 = rows2*cols2;
        
        % if this is only an update modify this parameters
        num_old_DAPI = 0;
        lfn = [filename.resultfolder filename.resultfile];
        load(lfn)
       
    
        Fname =  fieldnames(Results{1});
        Fname(1:2) = [];
        for thisiscoreclassification = 1
            try
                [~,~,CoreAssign] = xlsread([filename.datafolder 'corenumbercorrection.xlsx']);
            catch
                continue
            end
            CoreAssign =  cell2mat(CoreAssign);
            [~,~,TMAclasslist] = xlsread("Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\Results\HTMA399_Diagnosis_Map_DeID_v4.xls");
            oldresult = Results;
            [maxrow,~] = size(TMAclasslist);
            resultindexcol{1} = 'Resultindex';
            for i = 1:maxrow-1
                index = find(CoreAssign(:,2)==i);
                resultindexcol{i} = index;
            end
            TMAclasslist_new = [TMAclasslist(:,1:2),[{'ResultIndex'};resultindexcol.'],TMAclasslist(:,3:end)];
            for i = 3:2:maxrow
                TMAclasslist_new(i,4:33) = TMAclasslist_new(i-1,4:33);
            end
            clear i maxcol maxrow
            
            channel = 10;
            diagnosis = 'Glioblastoma';
            
            Results = oldresult;
            [maxrow,~] = size(TMAclasslist);
            markfordel = [];
            markforpos = [];
            poscount = 0;
            for i = 2:maxrow
                cdiag = TMAclasslist_new{i,channel}; 
                if isnan(cdiag)
                    cdiag = 'NAN';
                end
                recres = TMAclasslist_new{i,12}; 
                IDH = TMAclasslist_new{i,17}; 
                if ~contains(cdiag,diagnosis) || IDH ~= 0  || recres ~= 0 
                    resultindex = TMAclasslist_new{i,3};
                    TMAclasslist_new{i,34} = 1;
                    markfordel = [markfordel resultindex]; %markfordel marks all the non-IDH WT cores
                else
                    resultindex = TMAclasslist_new{i,3};
                    TMAclasslist_new{i,34} = 0;
                    markforpos = [markforpos resultindex];
                    poscount = [poscount i];
                end
            end
            Results(markfordel) = {''}; 
            clear i maxrow maxcol channel i diagnosis background cdiag channel cols2 cut i1 IDH index lfn maxrow maxtile2 num_old_DAPI
        end
        
        %roi selection
        try
            [~,~,CoreROI] = xlsread([filename.datafolder 'ROImask.xlsx']);
            CoreROI =  cell2mat(CoreROI);
            ROImask = CoreROI(:,2) == 1;
            Results(~ROImask) = {''}; 
        end
    
        AggrResults = [];
        for f = 1:length(Fname)
            AggrResults.(Fname{f}) = [];
            for i = 1:length(Results)
                if ~isempty(Results{i})
                    if strcmp(Fname{f},'Name')
                    AggrResults.(Fname{f}) = [AggrResults.(Fname{f});repmat(str2num(Results{i}.(Fname{f})),length(Results{i}.(Fname{f-1})),1)];
                    else
                      AggrResults.(Fname{f}) = [AggrResults.(Fname{f});Results{i}.(Fname{f})];
                    end
                end
            end
        end
        clear f i Fname
    
    
        nonimmune = AggrResults.MeanNucSign(:,1)>0;
        immunechannels = [52,54,40,46];
        immunecutoff = [100,100,1000,300];
        for i = 1:length(immunechannels)
            nonimmune = nonimmune & (AggrResults.MeanNucSign(:,immunechannels(i)) < immunecutoff(i));
        end
    
        Channelcutoffcell = readcell('Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Inteferons\ANALYSIS03302023\Channelscutoff.xlsx');
        [maxchannel,~] = size(Channelcutoffcell);
        AllChannel = Channelcutoffcell(:,1);
        AllChanneln = 1:maxchannel;
        Channelcutoff = cell2mat(Channelcutoffcell(:,2));
        amps = cell2mat(Channelcutoffcell(:,3));
        Channelchar = Channelcutoffcell(:,4);
        ChannelL = Channelcutoffcell(:,5);
        clear maxchannel Channelcutoffcell i1 lfn maxtile2 num_old_DAPI radius10X rows2 cols2 cut background foldinginitializationsteps
        
        % Neighborhood Analysis Prep
        AggrResults_micro = AggrResults;
        Results_micro = Results;
        load('Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Inteferons\ANALYSIS03302023\Results\2023-4-10_Results.mat') %Primary Nucleus Result
        Fname =  fieldnames(Results{1});
        Fname(1:2) = [];
        for f = 1:length(Fname)
            AggrResults.(Fname{f}) = [];
            for i = 1:length(Results)
                if ~isempty(Results{i})
                  AggrResults.(Fname{f}) = [AggrResults.(Fname{f});Results{i}.(Fname{f})];
                end
            end
        end
        AggrResults_primary = AggrResults;
        Results(markfordel) = {''}; 
        Results_primary = Results;
        AggrResults = AggrResults_micro;
        Results = Results_micro;
        
        Results_primarycorr = [];
         %Neighborhood analysis prep
        Results_primarycorr = [];
        cgaschannel = 42;
        for coren = 1:length(Results_micro)
            if isempty(Results_primary{coren})
                continue
            end
            if isempty(Results_micro{coren})
                Results_primary{coren}.micronucleus = Results_primary{coren}.Area<0;
                Results_primary{coren}.micronucleusrupture = Results_primary{coren}.Area<0;
                Results_primary{coren}.micronucleusmask_cgas = Results_primary{coren}.Area<0;
                continue
            end
            microCentroids = [Results_micro{coren}.CentroidX Results_micro{coren}.CentroidY];
            primaryCentroids = [Results_primary{coren}.CentroidX Results_primary{coren}.CentroidY];    
            micronucleusmask = primaryCentroids(:,1)<0;
            if ~isempty(primaryCentroids)
                [k,dist] = dsearchn(primaryCentroids,microCentroids);
                k(dist>100) = [];
                micronucleusmask(unique(k)) = true;
            end
            Results_primary{coren}.micronucleus = micronucleusmask;
    
            microrupturemask = (Results_micro{coren}.MeanNucSign(:,options.substructures.channel) > 300) & (Results_micro{coren}.MeanNucSign(:,options.substructures.channel) > 1.3*Results_micro{coren}.MeanCytSign(:,options.substructures.channel));
            microCentroids_rupture = [Results_micro{coren}.CentroidX(microrupturemask) Results_micro{coren}.CentroidY(microrupturemask)];
             micronucleusmask_rupture = primaryCentroids(:,1)<0;
            if ~isempty(microCentroids_rupture)
                [k_rupture,dist_rupture] = dsearchn(primaryCentroids,microCentroids_rupture);
                k_rupture(dist_rupture>100) = [];
                micronucleusmask_rupture(unique(k_rupture)) = true;
            end
            Results_primary{coren}.micronucleusrupture = micronucleusmask_rupture;
    
            cgasmask = (Results_micro{coren}.MeanNucSign(:,cgaschannel) > 200) & (Results_micro{coren}.MeanNucSign(:,cgaschannel) > 1.39*Results_micro{coren}.MeanCytSign(:,cgaschannel));
            microCentroids_cgas = [Results_micro{coren}.CentroidX(cgasmask) Results_micro{coren}.CentroidY(cgasmask)];
            micronucleusmask_cgas = primaryCentroids(:,1)<0;
            if ~isempty(microCentroids_cgas)
                [k_cgas,dist_cgas] = dsearchn(primaryCentroids,microCentroids_cgas);    
                k_cgas(dist_cgas>100) = [];
                micronucleusmask_cgas(unique(k_cgas)) = true;
            end
            Results_primary{coren}.micronucleusmask_cgas = micronucleusmask_cgas;
    
            Fname =  fieldnames(Results_primary{coren});
            for Fnamei = 1:length(Fname)
                orrg = Results_primary{coren}.(Fname{Fnamei});
                [mr,mc] = size(orrg);
                 if (mr > 1)
                        Results_primarycorr{coren}.(Fname{Fnamei}) = orrg(k,:);
                 elseif strcmp(Fname{Fnamei},'Name')
                        Results_primarycorr{coren}.(Fname{Fnamei}) = orrg;
                 else
                        Results_primarycorr{coren}.(Fname{Fnamei}) = orrg;
                 end
            end
        end
        Fname =  fieldnames(Results_primarycorr{end});
        AggrResults_primarycorr = [];
        for f = 1:length(Fname)
            AggrResults_primarycorr.(Fname{f}) = [];
            for i = 1:length(Results)
                if ~isempty(Results{i})
                  AggrResults_primarycorr.(Fname{f}) = [AggrResults_primarycorr.(Fname{f});Results_primarycorr{i}.(Fname{f})];
                end
            end
        end
        nonemptyindex = find(~cellfun(@isempty,Results_primary));
        Fname =  fieldnames(Results_primary{nonemptyindex(2)});
        AggrResults_primary = [];
        for f = 1:length(Fname)
            AggrResults_primary.(Fname{f}) = [];
            for i = 1:length(Results_primary)
                if ~isempty(Results_primary{i})
                  AggrResults_primary.(Fname{f}) = [AggrResults_primary.(Fname{f});Results_primary{i}.(Fname{f})];
                end
            end
        end
        processsave = [filename.resultfolder 'processsave414.mat'];
        save( processsave)
    end
    %%
    clear
    processsave = ['Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Inteferons\ANALYSIS03302023micro\Results\processsave414.mat'];
    load(processsave)
    %% NF kb test
    plist = [];
    count = 0;
    for coren = 1:length(Results)
        if isempty(Results{coren})
            plist(coren) = 9999;
            count = count+1;
            continue
        end
        cmeanBAFpos = Results{coren}.MeanNucSign(:,12);
        cmeanBAFneg = Results{coren}.MeanNucSign(:,12);
        cmeanBAFpos = cmeanBAFpos(Results{coren}.AreaSubstruct(:,1) > 0);
        cmeanBAFneg = cmeanBAFneg(Results{coren}.AreaSubstruct(:,1) == 0);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        plist(coren) = ci(1);
    end
    [~,index] = sort(plist);
    
    %% 
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & strcmp(Channelchar{i},'Interferon')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanBAFpos = AggrResults_primary.MeanFullCellSign(:,ccn);
        cmeanBAFneg = AggrResults_primary.MeanFullCellSign(:,ccn);
        cmeanBAFpos = cmeanBAFpos(AggrResults_primary.AreaSubstruct(:,1)>0 );
        cmeanBAFneg = cmeanBAFneg(AggrResults_primary.AreaSubstruct(:,1)==0 );
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults_primary.MeanNucSign(:,ccn);
        bafpos = AggrResults_primary.micronucleusrupture == true;
        bafneg = AggrResults_primary.micronucleusrupture == false;
        ksdensity(log10(aggrcurrent(bafpos)+1))
        ksdensity(log10(aggrcurrent(bafneg)+1))
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('BAF +','BAF -','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Micronucleus and Primary Nucleus')
        end
        hold off
    end
    % exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
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
    exportgraphics(gcf,'test.pdf','Resolution',600,'ContentType','vector')
    
    
    %% 
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & strcmp(Channelchar{i},'Interferon')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanBAFpos = AggrResults_primary.MeanFullCellSign(:,ccn);
        cmeanBAFneg = AggrResults_primary.MeanFullCellSign(:,ccn);
        cmeanBAFpos = cmeanBAFpos(AggrResults_primary.micronucleusrupture == true );
        cmeanBAFneg = cmeanBAFneg(AggrResults_primary.micronucleusrupture == false & AggrResults_primary.micronucleus == true );
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults_primary.MeanNucSign(:,ccn);
        bafpos = AggrResults_primary.micronucleusrupture == true;
        bafneg = AggrResults_primary.micronucleusrupture == false;
        ksdensity(log10(aggrcurrent(bafpos)+1))
        ksdensity(log10(aggrcurrent(bafneg)+1))
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('BAF +','BAF -','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Micronucleus and Primary Nucleus')
        end
        hold off
    end
    % exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
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
    exportgraphics(gcf,'test.pdf','Resolution',600,'ContentType','vector')
    
    %% Mean Nuc with micro rupture vs primary rupture
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & strcmp(Channelchar{i},'Interferon')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanBAFpos = AggrResults_primary.MeanNucSign(:,ccn);
        cmeanBAFneg = AggrResults_primary.MeanNucSign(:,ccn);
        cmeanBAFpos = cmeanBAFpos(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == true));
        cmeanBAFneg = cmeanBAFneg(AggrResults_primary.MeanFociSign(:,1)>0);
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults_primary.MeanNucSign(:,ccn);
        bafpos = AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == true);
        bafneg = AggrResults_primary.MeanFociSign(:,1)>0;
        ksdensity(log10(aggrcurrent(bafpos)+1))
        ksdensity(log10(aggrcurrent(bafneg)+1))
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('BAF +','BAF -','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Micronucleus and Primary Nucleus')
        end
        hold off
    end
    exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
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
    title('Difference in Primary Nucleus with Micronucleus Rutpure versus Primary Rupture')
    xlabel('Normalized Mean Nuclear Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    exportgraphics(gca,'test.jpg','Resolution',600)
    %% Mean Nuc with micro rupture vs no rupture, nuclear
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & strcmp(Channelchar{i},'Interferon') & strcmp(ChannelL{i},'N')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    
        
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        isNuc = ChannelL{channeln} == 'N';
        if isNuc
            soxchannel = 55;
            soxcutoff = 50;
            soxmask = (AggrResults_primary.MeanNucSign(:,soxchannel) > soxcutoff);
            cdata = AggrResults_primary.MeanNucSign(:,ccn);
            cmeanBAFpos = cdata(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == true) & soxmask);
            cmeanBAFneg = cdata(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == false) & soxmask);
    
        else
            cdata = AggrResults_micro.MeanCytSign(:,ccn);
            cgasmask = (AggrResults_micro.MeanNucSign(:,options.substructures.channel) > 300) & (AggrResults_micro.MeanNucSign(:,options.substructures.channel) > 1.3*AggrResults_micro.MeanCytSign(:,options.substructures.channel));
            cmeanBAFpos = cdata(cgasmask & soxmask);
            cmeanBAFneg = cdata(~cgasmask & soxmask);
        end
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults_primary.MeanNucSign(:,ccn);
        bafpos = AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == true);
        bafneg = AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == false);
        ksdensity(log10(aggrcurrent(bafpos)+1))
        ksdensity(log10(aggrcurrent(bafneg)+1))
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('BAF +','BAF -','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Micronucleus and Primary Nucleus')
        end
        hold off
    end
    % exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
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
            text(plotpoint,x(i),"**",'Color','black','FontSize',15)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',15)
        end
    end
    xlim([min(y)-0.3,max(y)+0.25])
    text(0.2,x(2),"*p < 0.05","FontSize",12)
    text(0.2,x(1),"**p < 0.01","FontSize",12)
    title('Difference in Primary Nucleus with Micronucleus Ruptured versus Unruptured')
    xlabel('Normalized Mean Nuclear Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    exportgraphics(gca,'test.jpg','Resolution',600)
    %% Mean Nuc with micro rupture vs no rupture, cyto
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & strcmp(Channelchar{i},'Interferon') & strcmp(ChannelL{i},'C')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
      
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        isNuc = ChannelL{channeln} == 'N';
        if isNuc
            cdata = AggrResults_primary.MeanNucSign(:,ccn);
            cmeanBAFpos = cdata(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == true) );
            cmeanBAFneg = cdata(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == false) );
        else
            soxchannel = 55;
            soxcutoff = 50;
            soxmask = (AggrResults_micro.MeanNucSign(:,soxchannel) > soxcutoff);
            cdata = AggrResults_micro.MeanCytSign(:,ccn);
            cgasmask = (AggrResults_micro.MeanNucSign(:,options.substructures.channel) > 300) & (AggrResults_micro.MeanNucSign(:,options.substructures.channel) > 1.3*AggrResults_micro.MeanCytSign(:,options.substructures.channel));
            cmeanBAFpos = cdata(cgasmask & soxmask);
            cmeanBAFneg = cdata(~cgasmask );
        end
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults_primary.MeanNucSign(:,ccn);
        bafpos = AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == true);
        bafneg = AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == false);
        ksdensity(log10(aggrcurrent(bafpos)+1))
        ksdensity(log10(aggrcurrent(bafneg)+1))
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('BAF +','BAF -','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Micronucleus and Primary Nucleus')
        end
        hold off
    end
    % exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
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
            text(plotpoint,x(i),"**",'Color','black','FontSize',15)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',15)
        end
    end
    xlim([min(y)-0.3,max(y)+0.25])
    text(0.2,x(2),"*p < 0.05","FontSize",12)
    text(0.2,x(1),"**p < 0.01","FontSize",12)
    title('Difference in Ruptured Micronucleus versus Unruptured')
    xlabel('Normalized Mean Cytoplasmic Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    exportgraphics(gca,'test.jpg','Resolution',600)
    %% Mean Nuc with micro cgas vs no cgas
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & strcmp(Channelchar{i},'Interferon')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        isNuc = ChannelL{channeln} == 'N';
        if isNuc
            cdata = AggrResults_primary.MeanNucSign(:,ccn);
            cmeanBAFpos = cdata(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == true));
            cmeanBAFneg = cdata(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusrupture == false));
        else
            cdata = AggrResults_micro.MeanCytSign(:,ccn);
            cgasmask = (AggrResults_micro.MeanNucSign(:,cgaschannel) > 200) & (AggrResults_micro.MeanNucSign(:,cgaschannel) > 1.39*AggrResults_micro.MeanCytSign(:,cgaschannel));
            cmeanBAFpos = cdata(cgasmask);
            cmeanBAFneg = cdata(~cgasmask);
        end
    
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    % exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
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
            text(plotpoint,x(i),"**",'Color','black','FontSize',15)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',15)
        end
    end
    xlim([min(y)-0.3,max(y)+0.25])
    text(0.2,x(2),"*p < 0.05","FontSize",12)
    text(0.2,x(1),"**p < 0.01","FontSize",12)
    title('Difference in Primary Nucleus with Micronucleus with CGAS versus wihtout CGAS')
    xlabel('Normalized Mean Nuclear Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    exportgraphics(gca,'test.jpg','Resolution',600)
    
    %% Mean Nuc with micro cgas+baf+ vs baf- cgas
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & strcmp(Channelchar{i},'Interferon')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanBAFpos = AggrResults_primary.MeanNucSign(:,ccn);
        cmeanBAFneg = AggrResults_primary.MeanNucSign(:,ccn);
        cmeanBAFpos = cmeanBAFpos(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusmask_cgas == true) & (AggrResults_primary.micronucleusrupture == true));
        cmeanBAFneg = cmeanBAFneg(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusmask_cgas == false) & (AggrResults_primary.micronucleusrupture == true));
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    % exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
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
            text(plotpoint,x(i),"**",'Color','black','FontSize',15)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',15)
        end
    end
    xlim([min(y)-0.3,max(y)+0.25])
    text(0.2,x(2),"*p < 0.05","FontSize",12)
    text(0.2,x(1),"**p < 0.01","FontSize",12)
    title('Difference in Primary Nucleus with Micronucleus with CGAS-BAF+ versus CGAS-BAF-')
    xlabel('Normalized Mean Nuclear Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    exportgraphics(gca,'test.jpg','Resolution',600)
    
    %% Mean Nuc with micro cgas+baf+ vs baf-cgas-
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & strcmp(ChannelL{i},'Interferon')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanBAFpos = AggrResults_primary.MeanNucSign(:,ccn);
        cmeanBAFneg = AggrResults_primary.MeanNucSign(:,ccn);
        cmeanBAFpos = cmeanBAFpos(AggrResults_primary.micronucleus == true & (AggrResults_primary.micronucleusmask_cgas == true) & (AggrResults_primary.micronucleusrupture == true));
        cmeanBAFneg = cmeanBAFneg(AggrResults_primary.micronucleus == true & ((AggrResults_primary.micronucleusmask_cgas == true) & (AggrResults_primary.micronucleusrupture == false)));
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    % exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
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
            text(plotpoint,x(i),"**",'Color','black','FontSize',15)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',15)
        end
    end
    xlim([min(y)-0.3,max(y)+0.25])
    text(0.2,x(2),"*p < 0.05","FontSize",12)
    text(0.2,x(1),"**p < 0.01","FontSize",12)
    title('Difference in Primary Nucleus with Micronucleus with CGAS+BAF+ versus CGAS+BAF-')
    xlabel('Normalized Mean Nuclear Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    exportgraphics(gca,'test.jpg','Resolution',600)
    %% Binary data
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~isempty(Channelchar{i}) & strcmp(Channelchar{i},'Interferon')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanmarkerpos = AggrResults_primary.MeanNucSign(:,ccn)>Channelcutoff(ccn);
        cbafpos = AggrResults_primary.MeanFociSign(:,ccn)>0;
        [tbl,chi2,p,labels] = crosstab(cmeanmarkerpos,cbafpos);
        resulttable{2,cn+1} = sum(cmeanmarkerpos&cbafpos)/sum(cmeanmarkerpos&~cbafpos)-sum(~cmeanmarkerpos&cbafpos)/sum(~cmeanmarkerpos&~cbafpos);
        resulttable{3,cn+1} = p;
    end
    % exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
    hold on
    y = [resulttable{2,2:end}];
    [y,index] = sort(y);
    xsorted = resulttable(1,2:end);
    xsorted = xsorted(index);
    x = categorical(xsorted);
    x = reordercats(x,xsorted);
    barh(x,y);
    %xlim([min(y)-0.3,max(y)+0.25])
    % text(0.2,x(2),"*p < 0.05","FontSize",12)
    % text(0.2,x(1),"**p < 0.01","FontSize",12)
    title('Difference in Primary Nucleus with Micronucleus Ruptured versus Unruptured')
    xlabel('Normalized Mean Nuclear Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    exportgraphics(gca,'test.jpg','Resolution',600)
    %% draw
    k1 =  268;
    fullstackfolder = [filename.datafolder 'dearray' filesep]; 
    FileTif = [fullstackfolder num2str(filename.realcoreinfo(k1).index) '.ome' filename.suffix];
    DAPI_img = uint16(imread(FileTif,'Index',1));
    imshow(DAPI_img*2)
    hold on
    for i = 1:length(microCentroids)
        plot([microCentroids(i,1),primaryCentroids(k(i),1)],[microCentroids(i,2),primaryCentroids(k(i),2)])
    end
    hold off
    
    % CorrPrimary Scimap prep
    c = [];
    for i = 1:length(Channelchar)
        try
            if (strcmp(Channelchar{i},'Nuclear') || strcmp(Channelchar{i},'Cytoplasmic')) && ChannelL{i} == 'L' %Palentirmarks{i} == 'Palantir' %
                c = [c i];
            end
        end
    end
    outL =(AggrResults_primarycorr.MeanNucSign(:,c));
    outL = log10(outL);
    [ar,ac]=find(~(outL>0));
    outL(ar,:) = [];
    
    
    T = array2table(outL);
    T.Properties.VariableNames(:) = AllChannel(c);
    writetable(T,[filename.resultfolder 'callprimcorr.csv'])
    
    focimask = AggrResults_primarycorr.MeanFociSign(:,1) > 0;
    meta = [];
    meta = array2table([AggrResults_primarycorr.Area,AggrResults_primarycorr.Solidity,AggrResults_primarycorr.CentroidX,AggrResults_primarycorr.CentroidY,focimask]);
    meta.Properties.VariableNames(:) = {'Area','Solidity','CentroidX','CentroidY','BAF'};
    meta(ar,:) = [];
    writetable(meta,[filename.resultfolder 'metaprimcorr.csv'])
    
    % CorrPrimary Scimap prep
    c = [];
    for i = 1:length(Channelchar)
        try
            if (strcmp(Channelchar{i},'Nuclear') || strcmp(Channelchar{i},'Cytoplasmic')) && ChannelL{i} == 'L' %Palentirmarks{i} == 'Palantir' %
                c = [c i];
            end
        end
    end
    outL =(AggrResults_primary.MeanFullCellSign(:,c));
    outL = log10(outL);
    [ar,ac]=find(~(outL>0));
    outL(ar,:) = [];
    gap = 1;
    outL = outL(1:gap:end,:);
    T = array2table(outL);
    T.Properties.VariableNames(:) = AllChannel(c);
    writetable(T,[filename.resultfolder 'primarynucleus.csv'])
    
    focimask = AggrResults_primary.MeanFociSign(:,1) > 0;
    meta = [];
    meta = array2table([AggrResults_primary.Corenum,AggrResults_primary.Area,AggrResults_primary.Solidity,AggrResults_primary.CentroidX,AggrResults_primary.CentroidY,focimask,AggrResults_primary.micronucleus]);
    meta.Properties.VariableNames(:) = {'Corenum','Area','Solidity','CentroidX','CentroidY','BAF','Micronucleus'};
    meta(ar,:) = [];
    meta = meta(1:gap:end,:);
    writetable(meta,[filename.resultfolder 'metaprimarynucleus.csv'])
    
    joined = [meta T];
    writetable(joined ,[filename.resultfolder 'AggregatedData_329_GBMTMA.csv'])
    
    %% Core level analysis
    AllCoresn = unique(AggrResults_primary.Corenum(:));
    Results_core = [];
    for Allcoresn_index = 1:length(AllCoresn)
        cCorenum = AllCoresn(Allcoresn_index);
        CurrentCoremask = AggrResults_primary.Corenum(:) == cCorenum;
        cResults = Results_primary{cCorenum};
        Results_core.MeanFullCellSign_average(Allcoresn_index,:) = mean(cResults.MeanFullCellSign,1);
        Results_core.MeanNucSign_average(Allcoresn_index,:) = mean(cResults.MeanNucSign,1);
        Results_core.MeanCytSign_average(Allcoresn_index,:) = mean(cResults.MeanCytSign,1);
        Results_core.BAFcounts(Allcoresn_index) = sum(cResults.MeanInsideSubstruct(:,1)>0); 
        Results_core.Cellcounts(Allcoresn_index) = sum(cResults.MeanInsideSubstruct(:,2)>0); 
        Results_core.micronucleus(Allcoresn_index) = sum(cResults.micronucleus);
        Results_core.micronucleusrupture(Allcoresn_index) = sum(cResults.micronucleusrupture);
        Results_core.micronucleusmask_cgas(Allcoresn_index) = sum(cResults.micronucleusmask_cgas);
    end
    Results_core.BAFrate = (Results_core.BAFcounts./Results_core.Cellcounts);
    c = [];
    for i = 1:length(Channelchar)
        if ~isempty(Channelchar{i}) & strcmp(Channelchar{i},'Interferon')
            c = [c i];
        end
    end
    close all
    CorN = ChannelL(c);
    for i = 1:length(c)
        figure
        currentc = c(i);
        x = Results_core.BAFrate;
        if CorN{i} == 'N'
            CorNtype = 'Nuclear';
            y = Results_core.MeanNucSign_average(:,currentc);
        elseif CorN{i} == 'C'
            CorNtype = 'Cytoplasmic';
            y = Results_core.MeanCytSign_average(:,currentc);
        end
        [x_sorted,index_sorted] = sort(x);
        y_sorted = y(index_sorted);
        plot(x_sorted,y_sorted)
        title([AllChannel{currentc} ''])
        xlabel('Primary Nucleus Rupture Rate')
        ylabel(['Average Mean ' CorNtype ' Intensity'])
        exportgraphics(gca,[num2str(i) '.jpg'],'Resolution',600)
    end
    %% 
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & strcmp(Channelchar{i},'BAFanditsfriends')
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanBAFpos = AggrResults_primary.MeanFullCellSign(:,ccn);
        cmeanBAFneg = AggrResults_primary.MeanFullCellSign(:,ccn);
        cmeanBAFpos = cmeanBAFpos(AggrResults_primary.AreaSubstruct(:,1)>0 );
        cmeanBAFneg = cmeanBAFneg(AggrResults_primary.AreaSubstruct(:,1)==0 );
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults_primary.MeanNucSign(:,ccn);
        bafpos = AggrResults_primary.micronucleusrupture == true;
        bafneg = AggrResults_primary.micronucleusrupture == false;
        ksdensity(log10(aggrcurrent(bafpos)+1))
        ksdensity(log10(aggrcurrent(bafneg)+1))
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('BAF +','BAF -','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Micronucleus and Primary Nucleus')
        end
        hold off
    end
    % exportgraphics(gcf,'screening.pdf','Resolution',600)
    clf
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
    exportgraphics(gcf,'test.pdf','Resolution',600,'ContentType','vector')
    
    %% Summary Mean Nuc
    clear resulttable
    resulttable{1,1}='Total';
    resulttable{1,2}= sum(AggrResults.AreaSubstruct(:,1) >= 0);
    resulttable{1,3}='BAFPositive';
    resulttable{1,4}= sum(AggrResults.AreaSubstruct(:,1) > 0);
    resulttable{3,1}='Total Positive in BAF Foci';
    resulttable{4,1}='Positive % in BAF Foci';
    resulttable{3,1}='Totalpositivewithbaf';
    resulttable{4,1}='Percentage with baf';
    resulttable{5,1}='Totalpositivewithoutbaf';
    resulttable{6,1}='Percentage without baf';
    c = [1 7 ];
    for cn = 1:length(c)
        i = c(cn);
        ccn = AllChanneln(c(cn));
        totalpositive = sum(AggrResults.MeanNucSign(:,ccn) > Channelcutoff(i) & (AggrResults.MeanInsideSubstruct(:,1)>0));
        resulttable{3,cn+1} = totalpositive;
        resulttable{4,cn+1} = totalpositive/sum(AggrResults.MeanInsideSubstruct(:,1)>0);
        totalnegative =  sum(AggrResults.MeanNucSign(:,ccn) > Channelcutoff(i) & (AggrResults.MeanInsideSubstruct(:,1)==0));
        resulttable{5,cn+1} = totalnegative;
        resulttable{6,cn+1} = totalnegative/sum(AggrResults.MeanInsideSubstruct(:,1)==0);
        resulttable{7,cn+1} = sum(AggrResults.MeanNucSign(:,ccn) <= Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)>0);
        resulttable{8,cn+1} = sum(AggrResults.MeanNucSign(:,ccn) <= Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)==0);
        %fishertest
        tempmatrix = [resulttable{3,cn+1},resulttable{7,cn+1};resulttable{5,cn+1},resulttable{8,cn+1}];
        [~,p]=fishertest(tempmatrix);
        resulttable{9,cn+1} = p;
    end
    
    clear ccn cn i c p totalnegative totalpositive tempmatrix 
end