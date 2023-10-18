function Dataset_2_Step5_ResulttoCylinter
    %Read result
    load([filename.resultfolder filename.resultfile])

    fullstackfolder = [filename.datafolder 'dearray' filesep]; 
    ilastikfolder = [filename.analfolder filename.ilastiksegfol];
    ilastiksegfolder = [ilastikfolder 'seg\'];
    ilastiklabelfolder = [ilastikfolder 'label\'];
    cyfol = [filename.analfolder filename.cylinterfol];
    cyinputfol = [filename.analfolder filename.cylinterfol filename.cylinterinputfol];
    cymaskfol = [filename.analfolder filename.cylinterfol filename.cylinterinputfol filename.cylintermaskfol];
    cysegfol = [filename.analfolder filename.cylinterfol filename.cylinterinputfol filename.cylintersegfol];
    cytiffol = [filename.analfolder filename.cylinterfol filename.cylinterinputfol filename.cylintertiffol];
    try
        mkdir(cyfol)
        mkdir(cyinputfol)
        mkdir(cymaskfol)
        mkdir(cysegfol)
        mkdir(cytiffol)
    end
    Fname =  fieldnames(Results{1});
    save([cyinputfol '\Results'],'Results')
    
    for thisiscoreclassification = 1
        [~,~,CoreAssign] = xlsread('Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Inteferons\corenumbercorrection.xlsx');
        CoreAssign =  cell2mat(CoreAssign);
        [~,~,TMAclasslist] = xlsread("Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\Results\HTMA399_Diagnosis_Map_DeID_v4.xls");
        oldresult = Results;
        [maxrow,~] = size(TMAclasslist);
        resultindexcol{1} = 'Resultindex';
        for i = 2:maxrow
            index = find(CoreAssign(:,2)==i-1);
            resultindexcol{i} = index;
        end
        TMAclasslist_new = [TMAclasslist(:,1:2),resultindexcol.',TMAclasslist(:,3:end)];
        for i = 3:2:maxrow
            TMAclasslist_new(i,4:33) = TMAclasslist_new(i-1,4:33);
        end
        clear i maxrow maxcol channel i diagnosis background cdiag cols2 cut i1 IDH index lfn maxrow maxtile2 num_old_DAPI
        for n = 0:length(Results)-1
            str = 'unmicst-number_MeanFullCellSign: ["number", "UnknownL", "UnknownS", "CANCER-TRUE", 1]';
            newstr = strrep(str,'number',num2str(n));
            for i = 2:length(resultindexcol)
                if resultindexcol{i} == n+1
                    currentrow = i;
                end
            end
            diag = TMAclasslist_new{currentrow,10};
            if ~isnan(diag)
                newstr = strrep(newstr,'UnknownL',diag);
                newstr = strrep(newstr,'UnknownS',diag(1:3));
            end
            fprintf(newstr)
            fprintf('\n')
        end
    end
    try
        [~,~,CoreROI] = xlsread([filename.datafolder 'ROImask.xlsx']);
        CoreROI =  cell2mat(CoreROI);
        ROImask = CoreROI(:,2) == 1;
        %Results(~ROImask) = {''}; 
        for n = 0:length(Results)-1
            if ROImask(n+1) == 1
                str = '"number",';
                newstr = strrep(str,'number',num2str(n));
                fprintf(newstr)
            end
        end
    end

    for i = 1:length(Results)
        ni = i -1;
        copyfile([fullstackfolder '\' Results{i}.Name '.ome.tif'],[cytiffol '\' num2str(ni) '.tif'])
    end   
    for i = 1:length(Results)
        ni = i -1;
        copyfile([ilastiklabelfolder '\' Results{i}.Name '_label.tif'],[cymaskfol '\' num2str(ni) '.tif'])
    end
    for i = 1:length(Results)
        ni = i -1;
        copyfile([ilastiksegfolder '\' Results{i}.Name '_Seg.tif'],[cysegfol '\' num2str(ni) '.tif'])
    end
end