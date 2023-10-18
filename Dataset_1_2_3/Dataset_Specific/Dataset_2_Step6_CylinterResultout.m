function Dataset_2_Step6_CylinterResultout(filename,options)
    load([filename.resultfolder filename.resultfile])
    cylinteroutfile =  xlsread('Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Inteferons\ANALYSIS03302023micro\Results\selectROIs_412.csv');
    %cylinteroutfile = readtable('Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2022_11_GBM399_TMA_Lineage\LSP14355_INCell\ANALYSIS20220218micro\Results\selectROIs_backup.csv');
    %Eliminate Nucleus
    Fname =  fieldnames(Results{1});
    for i = 1:length(Results)
        pythonindex = i-1;
        ccylintercells = find(cylinteroutfile(:,66) == pythonindex);
        ccellindex = cylinteroutfile(:,2);
        ccellindex = ccellindex(ccylintercells)+1;
        for fi = 1:length(Fname)
            fcname = Fname{fi};
            if contains(fcname,'Foci')
                continue
            end
            fcvar = Results{i}.(fcname);
            [mr,mc] = size(fcvar);
            if mr == 1
                newResults{i}.(fcname) = fcvar;
            elseif mc == 1
                newResults{i}.(fcname) = fcvar(ccellindex);
            else
                newResults{i}.(fcname) = fcvar(ccellindex,:);
            end
        end
    end
    
    for i = 1:length(newResults)
        ni = str2double(newResults{i}.Name);
        newnewResults{ni} = newResults{i};
    end
    Results = newnewResults;        

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

    save([filename.resultfolder '2023-4-10_Results_PostCylinter.mat'],'Results')
end