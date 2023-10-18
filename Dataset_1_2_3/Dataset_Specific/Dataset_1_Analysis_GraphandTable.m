function Dataset_1_Analysis_GraphandTable
    %% Initialization
    for foldinginitializationsteps=1
        clear 
        %Change All parameters below  
        filename.suffix = '.tif';
        filename.datafolder = 'Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\';
        filename.analfolder = [filename.datafolder 'ANALYSIS May 2022\'];
        filename.resultfolder = [filename.analfolder 'Results\'];
        filename.codesfolder = [filename.analfolder 'Codes Used\'];
        
        filename.wavelengths = {' wv UV - DAPI)',' wv Blue - FITC)',' wv Green - dsRed)',' wv Red - Cy5)'};
        filename.prefix1 = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O'}; %
        filename.prefix2 = linspace(1,25,25);
        filename.dim = ['%02d']; 
        filename.dim2 = ['%01d']; % this is for the number of tiles in a core
        filename.cycleprefix = 'Cycle';
        filename.realcorecycleprefix = 'Cycle1';
        filename.cycles = [1 2 3 4 5 6 7 8];  % cycles to analyse
        filename.maxcycle = 8;
        
        
        filename.ilastiksegfol = 'Ilastik Segmentation\';
        filename.ilastiksuffix = '_Probabilities.tif';
        filename.segsuffix = '_Seg.tif';
        
        
        options.nuc = 1;
        options.cyt = 2;
        options.backgr = 3;
        options.max_prob = 65535;
        options.bkgdprob_min = 0.2;
        options.cytprob_min  = 0.8;
        options.nucprob_min  = 0.33;
        options.pausetime = 0;
        
        
        count = 0;
        for i1 = 1:length(filename.prefix1)
            for i2 = 1:length(filename.prefix2)
                count = count + 1;
                filename.realcoreinfo{count}.name = [filename.prefix1{i1} ' - ' num2str(filename.prefix2(i2),filename.dim) ];
                filename.realcoreinfo{count}.tiles = [1 2 3 4];
            end
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
        options.cellsize = 35;
        options.cytosize = options.cellsize/5;
        options.corrthresh = 0.5;
        options.appendflag = 0;
        options.dates = '2023-1-19_';
        
        filename.resultfile = [options.dates 'Results.mat'];
        
        
        options.substructures.flag = 1;
        options.substructures.channel = 12;
        options.substructures.bkgsize = 15;
        options.substructures.SV_thresh = 3;
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
        
        clear background cols2 count cut DAPIslice i1 i2 maxtile2 num_old_DAPI radius10X rows2
    
        lfn = ['Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\Results\' filename.resultfile];
        load(lfn)
    
        Fname =  fieldnames(Results{1});
        for i = 1:length(Results)
            if ~isempty(Results{i})
                Results{i}.Corenum(1:size(Results{i}.Area,1),1) = i;
                if size(Results{i}.MeanFociSign,1) < size(Results{i}.Area,1) %Correction for NaN being deleted 
                    Results{i}.MedianFociSign( size(Results{i}.Area,1),1) = uint16(0);
                    Results{i}.MeanFociSign( size(Results{i}.Area,1),1) = 0;
                end
            end
        end
        Fname =  fieldnames(Results{1});
        for f = 1:length(Fname)
            AggrResults.(Fname{f}) = [];
            for i = 1:length(Results)
                if ~isempty(Results{i})
                  AggrResults.(Fname{f}) = [AggrResults.(Fname{f});Results{i}.(Fname{f})];
                end
            end
        end
        
        clear f i filename Fname
    end
    
    %% Neighborhood Analysis
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2','PU.1'};
    AllChanneln = [10,14,16,7,18,26,15,20,11,6];
    Channelcutoff = [0,0,0,150,350,300,160,50,120,300];
    amps = [0,0,0,1.3,1.3,0,0,0,0,1.3];
    
    
    c = [ 1 2 3];
    clf
    figure(1)
    sox2mask = AggrResults.MeanNucSign(:,11) > 120;
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults.MeanNucSign(:,ccn);
        bafpos = aggrcurrent(AggrResults.AreaSubstruct(:,1)>10 & sox2mask);
        bafneg = aggrcurrent(AggrResults.AreaSubstruct(:,1)==0 & sox2mask);
        ksdensity(log10(bafpos+1))
        ksdensity(log10(bafneg+1))
        [~,p] = ttest2(bafpos,bafneg);
        plist{i,1} = p;
        plist{i,2} = AllChannel{c(i)};
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('BAF +','BAF -','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Ruptured and Unruptured Nucleus')
        end
        hold off
        set(gcf,'Position',[400 400 500 300])
        xlim([0 6])
    end
    exportgraphics(gcf,'Lamins.pdf','Resolution',600,'ContentType','vector')
    clear aggrcurrent AllChannel AllChanneln amps bafneg bafpos c ccn Channelcutoff Channelnum i p
    %% Summary Mean Nuc
    clear resulttable
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2'};
    AllChanneln = [10,14,16,7,18,26,15,20,11];
    Channelcutoff = [0,0,0,150,350,300,160,50,120];
    amps = [0,0,0,1.3,1.3,0,0,0,0];
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
    c = [1 2 3 4 5 6 7 8 9];
    for cn = 1:length(c)
        i = c(cn);
        ccn = AllChanneln(c(cn));
        totalpositive = sum(AggrResults.MeanNucSign(:,ccn) > Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)>0);
        resulttable{3,cn+1} = totalpositive;
        resulttable{4,cn+1} = totalpositive/sum(AggrResults.MeanInsideSubstruct(:,1)>0);
        totalnegative =  sum(AggrResults.MeanNucSign(:,ccn) > Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)==0);
        resulttable{5,cn+1} = totalnegative;
        resulttable{6,cn+1} = totalnegative/sum(AggrResults.MeanInsideSubstruct(:,1)==0);
        resulttable{7,cn+1} = sum(AggrResults.MeanNucSign(:,ccn) <= Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)>0);
        resulttable{8,cn+1} = sum(AggrResults.MeanNucSign(:,ccn) <= Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)==0);
        %fishertest
        tempmatrix = [resulttable{3,cn+1},resulttable{7,cn+1};resulttable{5,cn+1},resulttable{8,cn+1}];
        [~,p]=fishertest(tempmatrix);
        resulttable{9,cn+1} = p;
    end
    
    clear ccn cn i c p totalnegative totalpositive tempmatrix AllChannel AllChanneln Channelcutoff Channelnum amps
    %% Summary Mean Nuc
    clear resulttable
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2','PU.1'};
    AllChanneln = [10,14,16,7,18,26,15,20,11,6];
    Channelcutoff = [0,0,0,150,350,300,160,50,120,300];
    amps = [0,0,0,1.3,1.3,0,0,0,0,1.3];
    
    
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
    c = [10];
    for cn = 1:length(c)
        i = c(cn);
        ccn = AllChanneln(c(cn));
        totalpositive = sum(AggrResults.MeanNucSign(:,ccn) > Channelcutoff(i) & AggrResults.MeanNucSign(:,ccn) > amps(i)*AggrResults.MeanCytSign(:,ccn) & AggrResults.MeanFociSign(:,ccn) > 0 );
        resulttable{3,cn+1} = totalpositive;
        resulttable{4,cn+1} = totalpositive/sum(AggrResults.MeanInsideSubstruct(:,1)>0);
        totalnegative =  sum(AggrResults.MeanNucSign(:,ccn) > Channelcutoff(i) & AggrResults.MeanNucSign(:,ccn) > amps(i)*AggrResults.MeanCytSign(:,ccn) & AggrResults.MeanFociSign(:,ccn) == 0 ); 
        resulttable{5,cn+1} = totalnegative;
        resulttable{6,cn+1} = totalnegative/sum(AggrResults.MeanInsideSubstruct(:,1)==0);
        resulttable{7,cn+1} = sum(AggrResults.MeanNucSign(:,ccn) <= Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)>0);
        resulttable{8,cn+1} = sum(AggrResults.MeanNucSign(:,ccn) <= Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)==0);
        %fishertest
        tempmatrix = [resulttable{3,cn+1},resulttable{7,cn+1};resulttable{5,cn+1},resulttable{8,cn+1}];
        [~,p]=fishertest(tempmatrix);
        resulttable{9,cn+1} = p;
    end
    
    clear ccn cn i c p totalnegative totalpositive tempmatrix AllChannel AllChanneln Channelcutoff Channelnum amps
    
    %% Summary Mean Foci
    for thisiscoreclassification = 1
    load('Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\Results\2022-5-26_Resultsonly.mat')
    [~,~,TMAclasslist] = xlsread("Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\Results\HTMA399_Diagnosis_Map_DeID_v4.xls");
    oldresult = Results;
    [maxrow,~] = size(TMAclasslist);
    resultindexcol{1} = 'Resultindex';
    for i = 2:maxrow
       resultindexcol{i} = rem((i-2)*25+1,375)+floor(((i-2)*25+1)/375);
       if i-1 > 360
           resultindexcol{i} = rem((i-2)*25+1,375)+floor(((i-2)*25+1)/375)+300;
       end
    end
    TMAclasslist_new = [TMAclasslist(:,1:2),resultindexcol.',TMAclasslist(:,3:end)];
    for i = 3:2:maxrow
        TMAclasslist_new(i,4:33) = TMAclasslist_new(i-1,4:33);
    end
    clear i maxcol maxrow
    
    channel = 10;
    diagnosis = 'Glioblastoma';
    
    Results = oldresult;
    [maxrow,~] = size(TMAclasslist);
    markfordel = [];
    for i = 2:maxrow
        cdiag = TMAclasslist_new(i,channel); 
        recres = TMAclasslist_new{i,12}; 
        IDH = TMAclasslist_new{i,17}; 
        if ~strcmp(diagnosis,cdiag)  || IDH ~= 0 || recres ~= 0
            resultindex = TMAclasslist_new{i,3};
            TMAclasslist_new{i,34} = 1;
            markfordel = [markfordel resultindex]; %markfordel marks all the non-IDH WT cores
        else
            TMAclasslist_new{i,34} = 0;
        end
    end
    Results(markfordel) = {''}; 
    Fname =  fieldnames(Results{1});
    for i = 1:length(Results)
        if ~isempty(Results{i})
            Results{i}.Corenum(1:size(Results{i}.Area,1),1) = i;
            if size(Results{i}.MeanFociSign,1) < size(Results{i}.Area,1) %Correction for NaN being deleted 
                Results{i}.MedianFociSign( size(Results{i}.Area,1),1) = uint16(0);
                Results{i}.MeanFociSign( size(Results{i}.Area,1),1) = 0;
            end
        end
    end
    Fname =  fieldnames(Results{1});
    for f = 1:length(Fname)
        AggrResults.(Fname{f}) = [];
        for i = 1:length(Results)
            if ~isempty(Results{i})
              AggrResults.(Fname{f}) = [AggrResults.(Fname{f});Results{i}.(Fname{f})];
            end
        end
    end
    clear i maxrow maxcol channel i diagnosis
    end
    
    clear resulttable
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2','PU.1'};
    AllChanneln = [10,14,16,7,18,26,15,20,11,6];
    Channelcutoff = [0,0,0,150,350,300,160,50,120,500];
    amps = [0,0,0,1.3,1.3,0,0,0,0,1.3];
    sox2mask = AggrResults.MeanNucSign(:,11) > 120;
    resulttable{1,1}='Total';
    resulttable{1,2}= sum(AggrResults.AreaSubstruct(:,1) >= 0 & sox2mask);
    resulttable{1,3}='BAFPositive';
    resulttable{1,4}= sum(AggrResults.AreaSubstruct(:,1) > 0 & sox2mask);
    resulttable{3,1}='Total Positive in BAF Foci';
    resulttable{4,1}='Positive % in BAF Foci';
    
    c = [4 5 6];
    for cn = 1:length(c)
        i = c(cn);
        ccn = AllChanneln(c(cn));
        amp = amps(i);
        all = (AggrResults.MeanFociSign(:,ccn) >0);
        resulttable{2,cn+1} =  AllChannel(c(cn));
        posmask = (AggrResults.MeanFociSign(:,ccn) > Channelcutoff(i)) & (AggrResults.MeanFociSign(:,ccn) > (AggrResults.MeanNucSign(:,ccn) * amp));
        resulttable{3,cn+1} = sum(posmask & sox2mask);
        resulttable{4,cn+1} = resulttable{3,cn+1}/sum(all & sox2mask);
    end
    
    resulttable = [];
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
    c = [10];
    for cn = 1:length(c)
        i = c(cn);
        ccn = AllChanneln(c(cn));
        resulttable{2,cn+1} =  AllChannel(c(cn));
        totalpositive = sum(AggrResults.MeanNucSign(:,ccn) > Channelcutoff(i) & AggrResults.MeanNucSign(:,ccn) > amps(cn)* AggrResults.MeanCytSign(:,ccn) & AggrResults.MeanInsideSubstruct(:,1)>0);
        resulttable{3,cn+1} = totalpositive;
        resulttable{4,cn+1} = totalpositive/sum(AggrResults.MeanInsideSubstruct(:,1)>0);
        totalnegative =  sum(AggrResults.MeanNucSign(:,ccn) > Channelcutoff(i) & AggrResults.MeanNucSign(:,ccn) > amps(cn)* AggrResults.MeanCytSign(:,ccn) & AggrResults.MeanInsideSubstruct(:,1)==0);
        resulttable{5,cn+1} = totalnegative;
        resulttable{6,cn+1} = totalnegative/sum(AggrResults.MeanInsideSubstruct(:,1)==0);
        resulttable{7,cn+1} = sum(AggrResults.MeanNucSign(:,ccn) <= Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)>0);
        resulttable{8,cn+1} = sum(AggrResults.MeanNucSign(:,ccn) <= Channelcutoff(i) & AggrResults.MeanInsideSubstruct(:,1)==0);
        %fishertest
        tempmatrix = [resulttable{3,cn+1},resulttable{7,cn+1};resulttable{5,cn+1},resulttable{8,cn+1}];
        [~,p]=fishertest(tempmatrix);
        resulttable{9,cn+1} = p;
    end
    
    clear AllChannel AllChanneln Channelcutoff Channelnum amps c ccn cn i amp
    
    %% Analysis by cores
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2','DAPI'};
    AllChanneln = [10,14,16,7,18,26,15,20,11,1];
    Channelcutoff = [0,0,0,150,200,300,160,50,120];
    amps = [0,0,0,1.3,1.3,0,0,0,0];
    Coi = [4 5];
    n = 0;
    ctotaltotal = 0;
    ctotalpos = 0;
    for i = 1:length(Results)
        for ci = 1:length(Coi)
            ccn = AllChanneln(Coi(ci));
            caggrmask = AggrResults.Corenum == i;
            cmeanfoci = AggrResults.MeanFociSign(:,ccn);
            cmeannuc = AggrResults.MeanNucSign(:,ccn);
            cmeanfoci = cmeanfoci(caggrmask);
            cmeannuc = cmeannuc(caggrmask);
            carea = AggrResults.AreaSubstruct(:,1);
            carea = carea(caggrmask);
            focimask = carea > 0;
            ctotal = sum(focimask);
            corestat(i,ci) = sum((cmeanfoci > Channelcutoff(Coi(ci))) & (cmeanfoci > (cmeannuc * amps(Coi(ci))))) / ctotal;
        end
        ctotalpos = ctotalpos + sum((cmeanfoci > Channelcutoff(Coi(ci))) & (cmeanfoci > (cmeannuc * amps(Coi(ci)))));
        ctotaltotal = ctotaltotal + ctotal;
        n = n+1; 
    end
    coinames = AllChannel(Coi);
    figure,
    title('CGAS and PH2AX at BAF Foci for Cores')
    ylabel('Mean Positive Percentage of BAF Foci(%)')
    boxplot(corestat*100,coinames)
    
    for thisiscoreclassification = 1
    load('Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\Results\2022-5-26_Resultsonly.mat')
    [~,~,TMAclasslist] = xlsread("Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\Results\HTMA399_Diagnosis_Map_DeID_v4.xls");
    oldresult = Results;
    [maxrow,~] = size(TMAclasslist);
    resultindexcol{1} = 'Resultindex';
    for i = 2:maxrow
       resultindexcol{i} = rem((i-2)*25+1,375)+floor(((i-2)*25+1)/375);
       if i-1 > 360
           resultindexcol{i} = rem((i-2)*25+1,375)+floor(((i-2)*25+1)/375)+300;
       end
    end
    TMAclasslist_new = [TMAclasslist(:,1:2),resultindexcol.',TMAclasslist(:,3:end)];
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
    for i = 2:maxrow
        cdiag = TMAclasslist_new(i,channel); 
        recres = TMAclasslist_new{i,12}; 
        IDH = TMAclasslist_new{i,17}; 
        if ~strcmp(diagnosis,cdiag) || IDH ~= 0 || recres ~= 0 
            resultindex = TMAclasslist_new{i,3};
            TMAclasslist_new{i,34} = 1;
            markfordel = [markfordel resultindex]; %markfordel marks all the non-IDH WT cores
        else
            markforpos = [markforpos resultindex];
            TMAclasslist_new{i,34} = 0;
        end
    end
    Results(markfordel) = {''}; 
    writecell(TMAclasslist_new,'temp.xlsx')
    clear i maxrow maxcol channel i diagnosis
    end
    
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2','DAPI'};
    AllChanneln = [10,14,16,7,18,26,15,20,11,1];
    Channelcutoff = [0,0,0,150,200,300,160,50,120];
    amps = [0,0,0,1.3,1.3,0,0,0,0];
    Coi = [4 5 6];
    n = 0;
    ctotaltotal = 0;
    ctotalpos = 0;
    for i = 1:length(Results)
        if ~any(i == markfordel)
            for ci = 1:length(Coi)
                ccn = AllChanneln(Coi(ci));
                caggrmask = AggrResults.Corenum == i;
                cmeanfoci = AggrResults.MeanFociSign(:,ccn);
                cmeannuc = AggrResults.MeanNucSign(:,ccn);
                cmeanfoci = cmeanfoci(caggrmask);
                cmeannuc = cmeannuc(caggrmask);
                carea = AggrResults.AreaSubstruct(:,1);
                carea = carea(caggrmask);
                focimask = carea > 0;
                ctotal = sum(focimask);
                corestat(i,ci) = sum((cmeanfoci > Channelcutoff(Coi(ci))) & (cmeanfoci > (cmeannuc * amps(Coi(ci))))) / ctotal;
            end
            ctotalpos = ctotalpos + sum((cmeanfoci > Channelcutoff(Coi(ci))) & (cmeanfoci > (cmeannuc * amps(Coi(ci)))));
            ctotaltotal = ctotaltotal + ctotal;
            n = n+1; 
        end
    end
    coinames = AllChannel(Coi);
    figure,
    boxplot(corestat*100,coinames)
    title('CGAS and PH2AX at BAF Foci for GBM Cores')
    ylabel('Mean Positive Percentage of BAF Foci(%)')
    exportgraphics(gcf,'CgasP2hAXcores.pdf','Resolution',600,'ContentType','vector')
    
    %% PU1
    for thisiscoreclassification = 1
    load('Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\Results\2022-5-26_Resultsonly.mat')
    [~,~,TMAclasslist] = xlsread("Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\Results\HTMA399_Diagnosis_Map_DeID_v4.xls");
    oldresult = Results;
    [maxrow,~] = size(TMAclasslist);
    resultindexcol{1} = 'Resultindex';
    for i = 2:maxrow
       resultindexcol{i} = rem((i-2)*25+1,375)+floor(((i-2)*25+1)/375);
       if i-1 > 360
           resultindexcol{i} = rem((i-2)*25+1,375)+floor(((i-2)*25+1)/375)+300;
       end
    end
    TMAclasslist_new = [TMAclasslist(:,1:2),resultindexcol.',TMAclasslist(:,3:end)];
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
    for i = 2:maxrow
        cdiag = TMAclasslist_new(i,channel); 
        recres = TMAclasslist_new{i,12}; 
        IDH = TMAclasslist_new{i,17}; 
        if ~strcmp(diagnosis,cdiag) || IDH ~= 0 || recres ~= 0 
            resultindex = TMAclasslist_new{i,3};
            TMAclasslist_new{i,34} = 1;
            markfordel = [markfordel resultindex]; %markfordel marks all the non-IDH WT cores
        else
            markforpos = [markforpos resultindex];
            TMAclasslist_new{i,34} = 0;
        end
    end
    Results(markfordel) = {''}; 
    writecell(TMAclasslist_new,'temp.xlsx')
    clear i maxrow maxcol channel i diagnosis
    end
    
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2','PU.1'};
    AllChanneln = [10,14,16,7,18,26,15,20,11,6];
    Channelcutoff = [0,0,0,150,350,300,160,50,120,500];
    amps = [0,0,0,1.3,1.3,0,0,0,0,1.3];
    c = [10];
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanBAFpos = AggrResults.MeanNucSign(:,ccn);
        cmeanBAFneg = AggrResults.MeanNucSign(:,ccn);
        mask = AggrResults.MeanFociSign(:,1) > 0;
        cmeanBAFpos = cmeanBAFpos(mask);
        cmeanBAFneg = cmeanBAFneg(~mask);
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p] = ttest2(cmeanBAFpos,cmeanBAFneg);
        resulttable{3,cn+1} = p;
    end
    
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults.MeanNucSign(:,ccn);
        bafpos = aggrcurrent(AggrResults.MeanInsideSubstruct(:,1)>0);
        bafneg = aggrcurrent(AggrResults.MeanInsideSubstruct(:,1)==0);
        ksdensity(log10(bafpos+1))
        ksdensity(log10(bafneg+1))
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('BAF +','BAF -','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Ruptured versus Unruptured Cells')
        end
        hold off
        set(gcf,'Position',[400 400 1000 300])
    end
    exportgraphics(gcf,'PU1-rupturevsnorupture.pdf','Resolution',600,'ContentType','vector')
    
    f2 = figure();
    %comparing lamin in pu+ with pu-
    c = [3];
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanpos = AggrResults.MeanNucSign(:,ccn);
        cmeanneg = AggrResults.MeanNucSign(:,ccn);
        posmask = AggrResults.MeanNucSign(:,6) > 500 & AggrResults.MeanNucSign(:,6) > 1.3* AggrResults.MeanCytSign(:,6);
        negmask = AggrResults.MeanNucSign(:,11) > 50 ;
        cmeanpos = cmeanpos(posmask);
        cmeanneg = cmeanneg(negmask);
        resulttable{2,cn+1} = mean(cmeanpos)/mean(cmeanneg);
        [~,p] = ttest2(cmeanBAFpos,cmeanneg);
        resulttable{3,cn+1} = p;
    end
    
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults.MeanNucSign(:,ccn);
        pos = aggrcurrent(posmask);
        neg = aggrcurrent(negmask);
        ksdensity(log10(pos+1))
        ksdensity(log10(neg+1))
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('PU.1+','Sox2+','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Myeloid and Tumor Cells')
        end
        hold off
        set(f2,'Position',[400 400 1000 300])
    end
    exportgraphics(gcf,'PU1-SOX2.pdf','Resolution',600,'ContentType','vector')
    %% Mean Nuc vs Mean Foci
    clear resulttable
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2','DAPI'};
    AllChanneln = [10,14,16,7,18,26,15,20,11,1];
    Channelcutoff = [0,0,0,150,350,300,160,50,120];
    c = [1:10];
    resulttable{2,1} = 'MeanFoci:MeanNuc';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        all = (AggrResults.MeanFociSign(:,ccn) >0);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanfocisign = AggrResults.MeanFociSign(:,ccn);
        cmeannucsign = AggrResults.MeanNucSign(:,ccn);
        mask = cmeanfocisign > 0 &  cmeannucsign > 0;
        cmeanfocisign = cmeanfocisign(mask);
        cmeannucsign = cmeannucsign(mask);
        resulttable{2,cn+1} = mean(cmeanfocisign./cmeannucsign);
        [~,p] = ttest2(cmeanfocisign,cmeannucsign);
        resulttable{3,cn+1} = p;
    end
    % for cn = 1:length(c)
    %     channeln = c(cn);
    %     ccn = AllChanneln(channeln);
    %     all = (AggrResults.MeanFociSign(:,ccn) >0);
    %     resulttable{1,cn+1} = AllChannel{channeln};
    %     cmeanfocisign = AggrResults.MeanFociSign(:,ccn);
    %     cmeanedgesign = AggrResults.MeanEdgeSign(:,ccn);
    %     mask = cmeanfocisign > 0 &  cmeanedgesign > 0;
    %     cmeanfocisign = cmeanfocisign(mask);
    %     cmeanedgesign = cmeanedgesign(mask);
    %     resulttable{2,cn+1} = mean(cmeanfocisign./cmeanedgesign);
    %     [~,p] = ttest2(cmeanfocisign,cmeanedgesign);
    %     resulttable{3,cn+1} = p;
    % end
    %% Mean Nuc vs BAF
    clear resulttable
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2','DAPI'};
    AllChanneln = [10,14,16,7,18,26,15,20,11,1];
    c = [3];
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        cmeanBAFpos = AggrResults.MeanNucSign(:,ccn);
        cmeanBAFneg = AggrResults.MeanNucSign(:,ccn);
        mask = AggrResults.MeanFociSign(:,1) > 0;
        cmeanBAFpos = cmeanBAFpos(mask);
        cmeanBAFneg = cmeanBAFneg(~mask);
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p] = ttest2(cmeanBAFpos,cmeanBAFneg);
        resulttable{3,cn+1} = p;
    end
    
    for i = 1:length(c)
        subplot(length(c),1,i)
        hold on
        ccn = AllChanneln(c(i));
        aggrcurrent = AggrResults.MeanNucSign(:,ccn);
        bafpos = aggrcurrent(AggrResults.MeanInsideSubstruct(:,1)>0);
        bafneg = aggrcurrent(AggrResults.MeanInsideSubstruct(:,1)==0);
        ksdensity(log10(bafpos+1))
        ksdensity(log10(bafneg+1))
        xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
        legend('BAF +','BAF -','Location','northeast')
        if i == 1
           title('Probability Density Distribution Comparison between Micronucleus and Primary Nucleus')
        end
        hold off
        set(gcf,'Position',[400 400 1000 300])
    end
    
    %% R2 for lamin a/c
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67','OLIG2','HIF1A','SOX2','DAPI'};
    AllChanneln = [10,14,16,7,18,26,15,20,11,1];
    
    % x = AggrResults.MeanNucSign(:,16);
    % y = AggrResults.MeanNucSign(:,12);
    % mask = x > 200 & y > 500; 
    % fitlm(x(mask),y(mask))
    
    %% Look into images
    DATAcells = [];
    AllChannel = {'LAMINB1','LAMINB2','LAMINAC','CGAS','PH2AX','KI67'};
    AllChanneln = [10,14,16,7,18,12,26];
    Channelcutoff = [0.4,0.25,0.55,0.5,0.7,0.5];
    amp = [0,0,0,1.3,0,0];
    i = 6;
    positives = table.Sample((table.(AllChannel{i}) > Channelcutoff(i)) & (table.cluster_2d(:) == 0) & areathreshold)
end