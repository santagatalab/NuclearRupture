function Dataset_3_Step6_CylinterResultout(filename,options)
    load([filename.resultfolder filename.resultfile])
    cylinteroutfile =  xlsread([filename.resultfolder 'selectROIs_510.csv']);
    %cylinteroutfile = readtable('Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2022_11_GBM399_TMA_Lineage\LSP14355_INCell\ANALYSIS20220218micro\Results\selectROIs_backup.csv');
    %Eliminate Nucleus
    Fname =  fieldnames(Results{end});
    for i = 1:length(Results)
        pythonindex = i-1;
        ccylintercells = find(cylinteroutfile(:,end-2) == pythonindex);
        ccellindex = cylinteroutfile(:,2);
        ccellindex = ccellindex(ccylintercells)+1;
        for fi = 1:length(Fname)
            fcname = Fname{fi};
            fcvar = Results{i}.(fcname);
            if isempty(fcvar)
                continue
            end
            [mr,mc] = size(fcvar);
            if mr == 1
                newResults{i}.(fcname) = fcvar;
            elseif mc == 1
                newResults{i}.(fcname) = fcvar(ccellindex);
            else
                newResults{i}.(fcname) = fcvar(ccellindex,:);
            end
        end
        if isempty(newResults{i}.Area)
            newResults{i} = [];
        end
    end
    Results = newResults;        

    Fname =  fieldnames(Results{1});
    for ri = 1:length(Results)
        for fi = 1:length(Fname)
            [mr,mc] = size(Results{1}.(Fname{fi}));
            if mr > 1 && mc > 2
                Cleanon = Fname{fi};
                mat = Results{ri}.(Cleanon);
                [mr,mc] = size(mat);
                for i = 1:mr
                    cycle = mc/4;
                    firstcycledapi = mat(i,1);
                    for c = 2:cycle
                        rdif = (mat(i,(c-1)*4+1)/(firstcycledapi));
                        if (rdif > 2) || (rdif < 0.3)
                            mat(i,(c-1)*4+2) = -1;
                            mat(i,(c-1)*4+3) = -1;
                            mat(i,(c-1)*4+4) = -1;
                        end
                    end
                end
                Results{ri}.(Cleanon) = mat;
            end
        end 
    end

    save([filename.resultfolder filename.resultfile(1:end-4) '_PostCylinter.mat'],'Results')
end