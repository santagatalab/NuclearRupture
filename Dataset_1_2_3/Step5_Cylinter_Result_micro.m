function Step5_Cylinter_Result_micro(filename,options)
    %Read result
    load([filename.resultfolder filename.resultfilemicro])

    % Create folder structure for cylinter
    fullstackfolder = filename.fullstackfolder;
    ilastikfolder = [filename.microNucfolder filename.ilastiksegfol];
    ilastiksegfolder = [ilastikfolder 'seg\'];
    ilastiklabelfolder = [ilastikfolder 'label\'];
    filename.analfolder = filename.microNucfolder;
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
    save([cyinputfol 'Results.mat'],'Results')

    % Generate individual file information
    for n = 0:length(Results)-1
        str = 'unmicst-number_MeanFullCellSign: ["number", "Glioblastoma", "GBM", "CANCER-TRUE", 1]';
        newstr = strrep(str,'number',num2str(n));
        fprintf(newstr)
        fprintf('\n')
    end

    % Transfer raw stacks to cylinter folder
    for i = 1:length(Results)
        ni = i -1;
        copyfile([fullstackfolder '\' Results{i}.Name ],[cytiffol '\' num2str(ni) '.tif'])
    end   
    % Transfer Labeled stacks to cylinter folder
    for i = 1:length(Results)
        ni = i -1;
        labelfile = [ilastiklabelfolder num2str(filename.realcoreinfo(i).index) '_label.tif'];
        copyfile(labelfile,[cymaskfol '\' num2str(ni) '.tif'])
    end
    % Transfer Segmented stacks to cylinter folder
    for i = 1:length(Results)
        ni = i -1;
        segfile = [ilastiksegfolder num2str(filename.realcoreinfo(i).index) filename.segsuffix];
        copyfile(segfile,[cysegfol '\' num2str(ni) '.tif'])
    end
end