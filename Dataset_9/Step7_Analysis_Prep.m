function Step7_Analysis_Prep(filename,options)
    %Aggregate DATA
    %Load Micronucleus data and aggregate them
    if exist([filename.resultfolder filename.resultfilemicro],"file")
        load([filename.resultfolder filename.resultfilemicro])
    else
        Results = [];
    end
    Fname =  fieldnames(Results{end});
    for f = 1:length(Fname)
        AggrResults.(Fname{f}) = [];
        AggrResults.corenum = [];
        for i = 1:length(Results)
            if ~isempty(Results{i}) && ~isempty(Results{i}.Area)
                if f == 1
                    AggrResults.corenum = [AggrResults.corenum;repmat(i,length(Results{i}.Area),1)];
                end
                [mr,mc] = size(Results{i}.(Fname{f}));
                if (mr > 1)
                    AggrResults.(Fname{f}) = [AggrResults.(Fname{f});Results{i}.(Fname{f})];
                end
            end
        end
    end
    Results_micro = Results;
    AggrResults_micro = AggrResults;
    clear Fname f AggrResults Results i mr mc    

    % Load Primary Nucleus data
    load([filename.resultfolder filename.resultfile])    
    Results_primary = Results;
    % Correlate micronucleus with the nearest primary nucleus
    Results_primarycorr = [];
    cgaschannel = options.cgaschannel;
    for coren = 1:length(Results_micro)
        % Check if Primary nucleus result is empty
        if isempty(Results_primary{coren})
            continue
        end
        % Check if there is any micronucleus 
        if isempty(Results_micro) || isempty(Results_micro{coren}) || isempty(Results_micro{coren}.Area) 
            Results_primary{coren}.micronucleus = Results_primary{coren}.Area<0;
            Results_primary{coren}.micronucleusrupture = Results_primary{coren}.Area<0;
            Results_primary{coren}.micronucleusmask_cgas = Results_primary{coren}.Area<0;
            continue
        end
        % Load Centroids
        microCentroids = [Results_micro{coren}.CentroidX Results_micro{coren}.CentroidY];
        primaryCentroids = [Results_primary{coren}.CentroidX Results_primary{coren}.CentroidY];    
        micronucleusmask = primaryCentroids(:,1)<0;
        % Use Nearest Point Search to correlate micronucleus to its primary
        % nucleus
        if ~isempty(primaryCentroids)
            [k,dist] = dsearchn(primaryCentroids,microCentroids);
            k(dist>100) = [];
            micronucleusmask(unique(k)) = true;
        end
        Results_primary{coren}.micronucleus = micronucleusmask;
        
        % Check if each micronucleus rupture
        microrupturemask = (Results_micro{coren}.MeanNucSign(:,options.substructures.channel) > 300) & (Results_micro{coren}.MeanNucSign(:,options.substructures.channel) > 1.3*Results_micro{coren}.MeanCytSign(:,options.substructures.channel));
        microCentroids_rupture = [Results_micro{coren}.CentroidX(microrupturemask) Results_micro{coren}.CentroidY(microrupturemask)];
         micronucleusmask_rupture = primaryCentroids(:,1)<0;
        if ~isempty(microCentroids_rupture)
            [k_rupture,dist_rupture] = dsearchn(primaryCentroids,microCentroids_rupture);
            k_rupture(dist_rupture>100) = [];
            micronucleusmask_rupture(unique(k_rupture)) = true;
        end
        Results_primary{coren}.micronucleusrupture = micronucleusmask_rupture;

        % Check if each micronucleus is a cgas positive
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

    % Aggreate Primary Nucleus Result
    Results = Results_primary;
    Fname =  fieldnames(Results{end});
    AggrResults.Corenum = [];
    for f = 1:length(Fname)
        AggrResults.(Fname{f}) = [];
        for i = 1:length(Results)
            if ~isempty(Results{i}) && ~isempty(Results{i}.Area)
                if f == 1
                    AggrResults.Corenum = [AggrResults.Corenum;repmat(i,length(Results{i}.Area),1)];
                end
                [mr,~] = size(Results{i}.(Fname{f}));
                if (mr > 1)
                    AggrResults.(Fname{f}) = [AggrResults.(Fname{f});Results{i}.(Fname{f})];
                end
            end
        end
    end
    try
        AggrResults = rmfield(AggrResults,'corenum');
    end
    try
        AggrResults = rmfield(AggrResults,'CoreFlag');
    end
    AggrResults_primary = AggrResults;
    clear f i
    % Save Aggregate MN and PN 
    save([filename.resultfolder 'Aggregate_' filename.resultfile],'AggrResults_micro','AggrResults_primary','filename','options','-v7.3') %Saving matrix

%     % Optional output csv file 
%     Channelcutoffcell = readcell([filename.datafolder 'Channelscutoff_BC.xlsx']);
%     [maxchannel,~] = size(Channelcutoffcell);
%     AllChannel = Channelcutoffcell(:,1);
%     AllChanneln = 1:maxchannel;
%     Channelcutoff = cell2mat(Channelcutoffcell(:,2));
%     amps = cell2mat(Channelcutoffcell(:,3));
%     Channelchar = Channelcutoffcell(:,4);
%     ChannelL = Channelcutoffcell(:,5);
%     ChannelC = Channelcutoffcell(:,6);
%     
%     % CorrPrimary Scimap prep
%     c = 1:length(AllChanneln);
%     outL =(AggrResults_primary.MeanFullCellSign(:,c));
%     T = array2table((1:length(outL))');
%     T.Properties.VariableNames(:) = {'CellID'};
%     T = [T array2table(outL)];
%     T.Properties.VariableNames(2:end) = AllChannel(c);
%     focimask = AggrResults_primary.MeanFociSign(:,1) > 0;
%     fociarea = log10(AggrResults_primary.AreaSubstruct(:,1));
%     fociarea(~(focimask>0)) =0;
%     meta = [];
%     meta = array2table([AggrResults_primary.Corenum,AggrResults_primary.Area,AggrResults_primary.Solidity,AggrResults_primary.FullCentroidX,AggrResults_primary.FullCentroidY,focimask,fociarea,AggrResults_primary.micronucleus,AggrResults_primary.micronucleusrupture,AggrResults_primary.micronucleusmask_cgas,AggrResults_primary.tumor]);
%     meta = [meta AggrResults_primary.type];
%     meta.Properties.VariableNames(:) = {'Corenum','Area','Solidity','X_centroid','Y_centroid','BAF','BAFarea','Micronucleus','Micronucleusrupture','Micronucleus_cgas','tumor','type'};
%     T = [T meta];
%     writetable(T,[filename.resultfolder 'Lineagev2_' filename.tissues{1} '.csv'])
end
