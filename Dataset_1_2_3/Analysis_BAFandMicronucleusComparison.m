function Analysis_BAFandMicronucleusComparison(filename,options)
    load([filename.resultfolder 'Aggregate_' filename.resultfile])

    Channelcutoffcell = readcell([filename.datafolder ]);
    [maxchannel,~] = size(Channelcutoffcell);
    AllChannel = Channelcutoffcell(:,1);
    AllChanneln = 1:maxchannel;
    Channelcutoff = cell2mat(Channelcutoffcell(:,2));
    amps = cell2mat(Channelcutoffcell(:,3));
    Channelchar = Channelcutoffcell(:,4);
    ChannelL = Channelcutoffcell(:,5);
    ChannelC = Channelcutoffcell(:,6);

    AggrResults = AggrResults_primary;


    %% Primary Nucleus with Rutpure versus no Rupture
    clear resulttable
    ChannelTypes = {'BAFfriends','Lineage','Immune'};
    for ct = 1:length(ChannelTypes)
        ChannelType = ChannelTypes{ct};
        c = [];
        for i = 1:length(Channelchar)
            if ~ismissing(Channelchar{i}) && strcmp(ChannelL{i},ChannelType)
                c = [c i];
            end
        end
        % c([2]) = [];
        resulttable = [];
        resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
        resulttable{3,1} = 't2test-pvalue';
        sox2mask = AggrResults_primary.MeanNucSign(:,15) > 50;
        tumormask = cellfun(@(c)strcmp(c,'Tumor'),AggrResults_primary.type);
        for cn = 1:length(c)
            channeln = c(cn);
            ccn = AllChanneln(channeln);
            resulttable{1,cn+1} = AllChannel{channeln};
            caggr = AggrResults_primary.MeanNucSign(:,ccn);%./AggrResults_primary.MeanNucSign(:,1);
            posmask = caggr>0;
            cposmask = AggrResults_primary.MeanInsideSubstruct(:,1)>10000;
            cnegmask = AggrResults_primary.MeanInsideSubstruct(:,1)==0;
            cmeanBAFpos = caggr(cposmask&sox2mask&posmask&tumormask);
            cmeanBAFneg = caggr(cnegmask&sox2mask&posmask&tumormask);
            resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
            [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
            error = ci/mean(cmeanBAFneg)+1;
            errlow(cn) = resulttable{2,cn+1} - error(1);
            errhigh(cn) = error(2) - resulttable{2,cn+1};
            resulttable{3,cn+1} = p;
        end
     
        figure,
        hold on
        y = [resulttable{2,2:end}]-1;
        [y,index] = sort(y);
        xsorted = resulttable(1,2:end);
        xsorted = xsorted(index);
        x = categorical(xsorted);
        x = reordercats(x,xsorted);
        errlow = errlow(index);
        barh(x,y);
        er = errorbar(y,x,errlow,'.',"horizontal");
        er.LineWidth = 1.5;
        er.Color = 'Blue';
        for i = 1:length(x)
            if resulttable{3,index(i)+1} < 0.01
                if y(i) > 0
                    plotpoint = y(i)+errlow(i)+0.01;
                else 
                    plotpoint = y(i)-errlow(i)-0.04;
                end
                text(plotpoint,x(i),"**",'Color','black','FontSize',12)
            elseif  resulttable{3,index(i)+1} < 0.05
                if y(i) > 0
                    plotpoint = y(i)+errlow(i)+0.01;
                else 
                    plotpoint = y(i)-errlow(i)-0.04;
                end
                text(plotpoint,x(i),'*','Color','black','FontSize',12)
            end
        end
        xlim([min(y)-0.2,max(y)+0.25])
        text(0.1,x(2),"*p < 0.05","FontSize",12)
        text(0.1,x(1),"**p < 0.01","FontSize",12)
        title(['Difference in Primary Nucleus with Rutpure versus no Rupture'])
        xlabel('Normalized Mean Nuclear Intensity Ratio')
        hold off
        set(gcf,'position',[400,-10,1000,600])
        exportgraphics(gca,['Difference in Primary Nucleus with Rutpure versus no Rupture in ', ChannelType,'.jpg'],'Resolution',600)
    end

   %% Primary Nucleus with Micronucleus versus with no Micronucleus
    clear resulttable
    ChannelTypes = {'BAFfriends','Lineage','Immune'};
    for ct = 1:length(ChannelTypes)
        ChannelType = ChannelTypes{ct};
        c = [];
        for i = 1:length(Channelchar)
            if ~ismissing(Channelchar{i}) && strcmp(ChannelL{i},ChannelType)
                c = [c i];
            end
        end
        resulttable = [];
        resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
        resulttable{3,1} = 't2test-pvalue';
        sox2mask = AggrResults_primary.MeanNucSign(:,15) > 50;
        for cn = 1:length(c)
            channeln = c(cn);
            ccn = AllChanneln(channeln);
            resulttable{1,cn+1} = AllChannel{channeln};
            caggr = AggrResults_primary.MeanNucSign(:,ccn);%./AggrResults_primary.MeanNucSign(:,1);
            posmask = caggr>0;
            cposmask = AggrResults_primary.micronucleus(:,1)>0;
            cnegmask = AggrResults_primary.micronucleus(:,1)==0;
            cmeanBAFpos = caggr(cposmask&sox2mask&posmask);
            cmeanBAFneg = caggr(cnegmask&sox2mask&posmask);
            resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
            [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
            error = ci/mean(cmeanBAFneg)+1;
            errlow(cn) = resulttable{2,cn+1} - error(1);
            errhigh(cn) = error(2) - resulttable{2,cn+1};
            resulttable{3,cn+1} = p;
        end
        
        figure,
        hold on
        y = [resulttable{2,2:end}]-1;
        [y,index] = sort(y);
        xsorted = resulttable(1,2:end);
        xsorted = xsorted(index);
        x = categorical(xsorted);
        x = reordercats(x,xsorted);
        errlow = errlow(index);
        barh(x,y);
        er = errorbar(y,x,errlow,'.',"horizontal");
        er.LineWidth = 1.5;
        er.Color = 'Blue';
        for i = 1:length(x)
            if resulttable{3,index(i)+1} < 0.01
                if y(i) > 0
                    plotpoint = y(i)+errlow(i)+0.01;
                else 
                    plotpoint = y(i)-errlow(i)-0.04;
                end
                text(plotpoint,x(i),"**",'Color','black','FontSize',12)
            elseif  resulttable{3,index(i)+1} < 0.05
                if y(i) > 0
                    plotpoint = y(i)+errlow(i)+0.01;
                else 
                    plotpoint = y(i)-errlow(i)-0.04;
                end
                text(plotpoint,x(i),'*','Color','black','FontSize',12)
            end
        end
        xlim([min(y)-0.2,max(y)+0.25])
        text(0.1,x(2),"*p < 0.05","FontSize",12)
        text(0.1,x(1),"**p < 0.01","FontSize",12)
        title(['Difference in Primary Nucleus with Micronucleus versus with no Micronucleus'])
        xlabel('Normalized Mean Nuclear Intensity Ratio')
        hold off
        set(gcf,'position',[400,-10,1000,600])
        exportgraphics(gca,['Difference in Primary Nucleus with Micronucleus versus with no Micronucleus in ', ChannelType,'.jpg'],'Resolution',600)
    end
    %% Primary Nucleus with Micronucleus rupture versus no Micronucleus rupture
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & Channelcutoff(i) > 0 & (strcmp(ChannelL(i),'Lineage') ) %InterferonLineage
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    sox2mask = AggrResults_primary.MeanNucSign(:,15) > 80;
    tumormask = AggrResults_primary.tumor == 1;
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
        aggrcurrent = AggrResults_primary.MeanFullCellSign(:,ccn);
        cposmask = AggrResults_primary.micronucleusrupture==1;
        cposmask_cores = unique(AggrResults_primary.Corenum(cposmask));
        cposmask_cores_mask = ismember(AggrResults_primary.Corenum,cposmask_cores); %cposmask_cores
        cnegmask = ~cposmask & AggrResults_primary.micronucleus==1;%& cposmask_cores_mask;
        cmeanBAFpos = aggrcurrent(cposmask&sox2mask&tumormask);
        cmeanBAFneg = aggrcurrent(cnegmask&sox2mask&tumormask);
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    figure,
    hold on
    y = [resulttable{2,2:end}]-1;
    [y,index] = sort(y);
    xsorted = resulttable(1,2:end);
    xsorted = xsorted(index);
    x = categorical(xsorted);
    x = reordercats(x,xsorted);
    errlow = errlow(index);
    barh(x,y);
    er = errorbar(y,x,errlow,'.',"horizontal");
    er.LineWidth = 1.5;
    er.Color = 'Blue';
    for i = 1:length(x)
        if resulttable{3,index(i)+1} < 0.01
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),"**",'Color','black','FontSize',12)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',12)
        end
    end
    xlim([min(y)-0.2,max(y)+0.25])
    text(0.1,x(2),"*p < 0.05","FontSize",12)
    text(0.1,x(1),"**p < 0.01","FontSize",12)
    title('Difference in Primary Nucleus with Micronucleus rupture versus no Micronucleus rupture')
    xlabel('Normalized Mean Nuclear and Cytoplasmic Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    %exportgraphics(gca,'test.jpg','Resolution',600)
    exportgraphics(gcf,'Difference in Primary Nucleus with Micronucleus rupture versus no Micronucleus rupture.pdf','Resolution',600,'ContentType','vector')

    %% Primary Nucleus with Rupture versus no Rupture
    clear resulttable
    c = [];
    for i = 1:length(Channelchar)
        if ~ismissing(Channelchar{i}) & Channelcutoff(i) > 0 & (strcmp(ChannelL(i),'Lineage') ) %InterferonLineage
            c = [c i];
        end
    end
    resulttable{2,1} = 'MeanBAFpos:MeanBAFneg';
    resulttable{3,1} = 't2test-pvalue';
    sox2mask = AggrResults_primary.MeanNucSign(:,15) > 80;
    tumormask = AggrResults_primary.tumor == 1;
    for cn = 1:length(c)
        channeln = c(cn);
        ccn = AllChanneln(channeln);
        resulttable{1,cn+1} = AllChannel{channeln};
    %     if Channelchar{channeln} == 'C'
    %       aggrcurrent = AggrResults_primary.MeanCytSign(:,ccn);
    %     else
    %       aggrcurrent = AggrResults_primary.MeanNucSign(:,ccn);
    %     end
        aggrcurrent = AggrResults_primary.MeanFullCellSign(:,ccn);
        cposmask = AggrResults_primary.AreaSubstruct(:,1) > 0;
        cnegmask = ~cposmask ;%& cposmask_cores_mask;
        posmask = aggrcurrent > 0;
        cmeanBAFpos = aggrcurrent(cposmask&sox2mask&tumormask&posmask);
        cmeanBAFneg = aggrcurrent(cnegmask&sox2mask&tumormask&posmask);
        resulttable{2,cn+1} = mean(cmeanBAFpos)/mean(cmeanBAFneg);
        [~,p,ci,stats] = ttest2(cmeanBAFpos,cmeanBAFneg);
        error = ci/mean(cmeanBAFneg)+1;
        errlow(cn) = resulttable{2,cn+1} - error(1);
        errhigh(cn) = error(2) - resulttable{2,cn+1};
        resulttable{3,cn+1} = p;
    end
    
    figure,
    hold on
    y = [resulttable{2,2:end}]-1;
    [y,index] = sort(y);
    xsorted = resulttable(1,2:end);
    xsorted = xsorted(index);
    x = categorical(xsorted);
    x = reordercats(x,xsorted);
    errlow = errlow(index);
    barh(x,y);
    er = errorbar(y,x,errlow,'.',"horizontal");
    er.LineWidth = 1.5;
    er.Color = 'Blue';
    for i = 1:length(x)
        if resulttable{3,index(i)+1} < 0.01
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),"**",'Color','black','FontSize',12)
        elseif  resulttable{3,index(i)+1} < 0.05
            if y(i) > 0
                plotpoint = y(i)+errlow(i)+0.01;
            else 
                plotpoint = y(i)-errlow(i)-0.04;
            end
            text(plotpoint,x(i),'*','Color','black','FontSize',12)
        end
    end
    xlim([min(y)-0.2,max(y)+0.25])
    text(0.1,x(2),"*p < 0.05","FontSize",12)
    text(0.1,x(1),"**p < 0.01","FontSize",12)
    title('Difference in Primary Nucleus with Rupture versus no Rupture')
    xlabel('Normalized Mean Nuclear and Cytoplasmic Intensity Ratio')
    hold off
    set(gcf,'position',[400,-10,1000,600])
    %exportgraphics(gca,'test.jpg','Resolution',600)
    exportgraphics(gcf,'Difference in Primary Nucleus with Rupture versus no Rupture.pdf','Resolution',600,'ContentType','vector')
end

