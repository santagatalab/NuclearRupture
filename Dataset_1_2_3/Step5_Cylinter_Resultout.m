function Step5_Cylinter_Resultout(filename,options)
    % Load Result File from Step 4
    load([filename.resultfolder filename.resultfile])
    % Load Cylinter output file
    cylinteroutfile =  xlsread([filename.resultfolder 'selectROIs_53.csv']);

    %Eliminate Nucleus not included in Cylinter output files
    Fname =  fieldnames(Results{1});
    for i = 1:length(Results)
        pythonindex = i-1;
        ccylintercells = cylinteroutfile(:,end-2) == pythonindex;
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
    
    % Fix Name
    for i = 1:length(newResults)
        ni = str2double(newResults{i}.Name);
        newnewResults{ni} = newResults{i};
    end
    Results = newnewResults;        

    % Optional check for cross-cycle registration
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
                        % Remove any cells with intensity difference greater 
                        % than 2*Firstcycledapi or less than 0.3*Firstcycledapi
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

    %
    save([filename.resultfolder 'PostCylinter_' filename.resultfile],'Results')
end