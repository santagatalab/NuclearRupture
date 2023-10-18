function Dataset_3_Analysis_LineageComparison(filename,options)
    %% Load Data
    % Load Result file
    lfn = [filename.resultfolder filename.resultfile];
    load(lfn)

    % Read Core assignment file
    [~,~,CoreAssign] = xlsread([filename.analfolder filename.CoreAssign]);
    CoreAssign =  cell2mat(CoreAssign);

    % Read TMA file information
    [~,~,TMAclasslist] = xlsread([filename.analfolder filename.TMAclasslist]);
    oldresult = Results;
    [maxrow,~] = size(TMAclasslist);
    resultindexcol{1} = 'Resultindex';
    for i = 2:maxrow
        index = find(CoreAssign(:,2)==i);
        resultindexcol{i} = index;
    end
    TMAclasslist_new = [TMAclasslist(:,1:2),resultindexcol.',TMAclasslist(:,3:end)];
    for i = 3:2:maxrow
        TMAclasslist_new(i,4:33) = TMAclasslist_new(i-1,4:33);
    end
    clear i maxcol maxrow
    
    % Channel number for diagnosis information
    channel = 10;

    % Diagnosis Information
    diagnosis_wanted = 'Glioblastoma';
    % IDH+ (1) or IDH- (0) 
    IDH_wanted = 0; 
    % Recurrent (1) or Primary (0)
    recres_wanted = 0;
    
    % Identify Cores of Interest
    Results = oldresult;
    [maxrow,~] = size(TMAclasslist);
    markfordel = [];
    markforpos = [];
    for i = 2:maxrow
        cdiag = TMAclasslist_new(i,channel); 
        recres = TMAclasslist_new{i,12}; 
        IDH = TMAclasslist_new{i,17}; 
        resultindex = TMAclasslist_new{i,3};
        if contains(diagnosis_wanted,cdiag)  && any(IDH == IDH_wanted) && any(recres == recres_wanted) 
            TMAclasslist_new{i,34} = 0;
            markforpos = [markforpos resultindex];
        else
            TMAclasslist_new{i,34} = 1;
            markfordel = [markfordel resultindex]; 
        end
    end
    Results(markfordel) = {''}; 
    clear i maxrow maxcol channel i diagnosis_wanted background cdiag channel cols2 cut i1 IDH index lfn maxrow maxtile2 num_old_DAPI
    
    % Create Aggregate Results
    clear AggrResults
      Fname =  fieldnames(Results{2});
      for f = 1:length(Fname)
        AggrResults.(Fname{f}) = [];
        for i = 1:length(Results)
            if ~isempty(Results{i})
               [mr,mc] = size(Results{i}.(Fname{f}));
               if (mr > 1)
                    AggrResults.(Fname{f}) = [AggrResults.(Fname{f});Results{i}.(Fname{f})];
               elseif strcmp(Fname{f},'Name')
                    AggrResults.(Fname{f}) = [AggrResults.(Fname{f});repmat(str2double(Results{i}.(Fname{f})),length(Results{i}.Area),1)];
               else
                    AggrResults.(Fname{f}) = [AggrResults.(Fname{f});repmat(Results{i}.(Fname{f}),length(Results{i}.Area),1)];
               end
            end
        end
      end
    clear f i Fname

    % Load Channel names, number, and cutoff Information
    Channelcutoffcell = readcell([filename.analfolder filename.Channelcutoff]);
    [maxchannel,~] = size(Channelcutoffcell);
    AllChannel = Channelcutoffcell(:,1);
    AllChanneln = 1:maxchannel;
    Channelcutoff = cell2mat(Channelcutoffcell(:,2));
    amps = cell2mat(Channelcutoffcell(:,3));
    Channelchar = Channelcutoffcell(:,4);
    ChannelL = Channelcutoffcell(:,5);
    ChannelC = Channelcutoffcell(:,6);
    clear maxchannel Channelcutoffcell i1 lfn maxtile2 num_old_DAPI radius10X rows2 cols2 cut background foldinginitializationsteps i mc mr recres resultindex dapislice poscount 

    % Neighborhood Analysis
    AggrResults_micro = AggrResults;
    Results_micro = Results;
    load([filename.primaryNucfolder 'Results\20230703_Results_PostCylinter.mat']) %Primary Nucleus Result
    Results(markfordel) = {''};
    Results_primary = Results;
    AggrResults = AggrResults_micro;
    Results = Results_micro;
    
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
    end

    Fname =  fieldnames(Results_primary{2});
    AggrResults_primary = [];
    for f = 1:length(Fname)
        AggrResults_primary.(Fname{f}) = [];
        for i = 1:length(Results_primary)
            if any(markforpos==i)
                if ~isempty(Results_primary{i})
                    AggrResults_primary.(Fname{f}) = [AggrResults_primary.(Fname{f});Results_primary{i}.(Fname{f})];
                end
            end
        end
    end
