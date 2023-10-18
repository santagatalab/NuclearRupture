function Dataset_1_Step5_Analysis_Prep(filename,options)
    %Aggregate DATA
    
    %Load Micronucleus data and aggregate them
    load(['Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSISmicroaug2022\Results\2022-10-7_Results.mat'])
    Fname =  fieldnames(Results{end});
    for f = 1:length(Fname)
        AggrResults.(Fname{f}) = [];
        AggrResults.corenum = [];
        for i = 1:length(Results)
            if ~isempty(Results{i})
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

    Results_primary = Results;
    % Correlate micronucleus with the nearest primary nucleus
    Results_primarycorr = [];
    cgaschannel = 7;
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

        microrupturemask = (Results_micro{coren}.MeanNucSign(:,options.substructures.channel) > 200) & (Results_micro{coren}.MeanNucSign(:,options.substructures.channel) > 1.2*Results_micro{coren}.MeanCytSign(:,options.substructures.channel));
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

    Results = Results_primary;
    Fname =  fieldnames(Results{end});
    for f = 1:length(Fname)
        AggrResults.(Fname{f}) = [];
        AggrResults.corenum = [];
        for i = 1:length(Results)
            if ~isempty(Results{i})
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
    AggrResults_primary = AggrResults;
    clear f i

    save([filename.resultfolder 'Aggregate_' filename.resultfile],'AggrResults_micro','AggrResults_primary','filename','options','-v7.3') %Saving matrix

    %     Fname =  fieldnames(Results{end});
    %     for f = 1:length(Fname)
    %         AggrResults.(Fname{f}) = [];
    %         for i = 1:length(Results)
    %             if f == 1
    %                 AggrResults.corenum = [AggrResults.(Fname{f});repmat(i,length(Results{i}.Area),1)];
    %             end
    %             if ~isempty(Results{i})
    %                 [mr,mc] = size(Results{i}.(Fname{f}));
    %                 if (mr > 1)
    %                     AggrResults.(Fname{f}) = [AggrResults.(Fname{f});Results{i}.(Fname{f})];
    %                 end
    %             end
    %         end
    %     end
    %     Results_primary = Results;
    %     AggrResults_primary = AggrResults;
    %     clear Fname f AggrResults Results i mr mc 

            %Test micronucleus neighboor 
    %     k1 = 219;
    %     fullstackfolder = [filename.datafolder 'dearray' filesep]; 
    %     FileTif = [fullstackfolder num2str(filename.realcoreinfo(k1).index) '.ome' filename.suffix];
    %     DAPI_img = uint16(imread(FileTif,'Index',1));
    %     coren = k1;
    %     microCentroids = [Results_micro{coren}.CentroidX Results_micro{coren}.CentroidY];
    %     primaryCentroids = [Results_primary{coren}.CentroidX Results_primary{coren}.CentroidY];    
    %     [k,dist] = dsearchn(primaryCentroids,microCentroids);
    %     distupperlimit = 80;
    %     k(dist<distupperlimit) = [];
    %     imshow(DAPI_img*10)
    %     hold on
    %     microCentroids(dist<distupperlimit,:) = []; 
    %     for i = 1:length(microCentroids)
    %         plot([microCentroids(i,1),primaryCentroids(k(i),1)],[microCentroids(i,2),primaryCentroids(k(i),2)])
    %     end
    %     hold off
end
