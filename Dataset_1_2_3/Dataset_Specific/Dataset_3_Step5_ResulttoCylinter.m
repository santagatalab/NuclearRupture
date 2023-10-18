function Dataset_3_Step5_ResulttoCylinter(filename)
%Read result
    addpath(filename.codesfolder)
    load([filename.resultfolder filename.resultfile])

    fsfol = [filename.analfolder filename.fullstackfol];
    ilfol = [filename.analfolder filename.ilastiksegfol];
    illabfol = [filename.analfolder 'label\'];
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
    i = 1;
    oi = 1;
    while i <= length(Results)
        if ~isempty(Results{i}.Area)
            Results{i}.Corenum(1:size(Results{i}.Area,1),1) = oi;
            if size(Results{i}.MeanFociSign,1) < size(Results{i}.Area,1) %Correction for NaN being deleted
                Results{i}.MedianFociSign( size(Results{i}.Area,1),1) = uint16(0);
                Results{i}.MeanFociSign( size(Results{i}.Area,1),1) = 0;
            end
            i = i+1;
            oi = oi+1;
        else
            Results(i) = [];
            oi = oi+1;
        end
    end
    save([cyinputfol '\Results'],'Results')

    for n = 0:length(Results)-1
        str = 'unmicst-number_MeanNucSign: ["number", "Glioblastoma", "GBM", "CANCER-TRUE", 1]';
        newstr = strrep(str,'number',num2str(n));
        fprintf(newstr)
        fprintf('\n')
    end

    for n = 0:length(Results)-1
        str = '"number",';
        newstr = strrep(str,'number',num2str(n));
        fprintf(newstr)
    end

    for i = 1:length(Results)
        ni = i -1;
        copyfile([fsfol '\' Results{i}.Name '.ome.tif'],[cytiffol '\' num2str(ni) '.tif'])
    end
    for i = 1:length(Results)
        ni = i -1;
        copyfile([illabfol '\' Results{i}.Name '_label.tif'],[cymaskfol '\' num2str(ni) '.tif'])
    end
    for i = 1:length(Results)
        ni = i -1;
        copyfile([ilfol '\' Results{i}.Name '_Seg.tif'],[cysegfol '\' num2str(ni) '.tif'])
    end
end
