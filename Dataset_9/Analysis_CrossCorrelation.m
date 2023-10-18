function [outputArg1,outputArg2] = Analysis_CrossCorrelation(inputArg1,inputArg2)
    load([filename.resultfolder 'Aggregate_' filename.resultfile])

    Channelcutoffcell = readcell([filename.datafolder 'Channelscutoff_BC.xlsx']);
    [maxchannel,~] = size(Channelcutoffcell);
    AllChannel = Channelcutoffcell(:,1);
    AllChanneln = 1:maxchannel;
    Channelcutoff = cell2mat(Channelcutoffcell(:,2));
    amps = cell2mat(Channelcutoffcell(:,3));
    Channelchar = Channelcutoffcell(:,4);
    ChannelL = Channelcutoffcell(:,5);
    ChannelC = Channelcutoffcell(:,6);

    AggrResults = AggrResults_primary;
    %% Cross Correlation
    c = [];
    for i = 1:length(Channelchar)
        try
            if ~ismissing(ChannelL{i}) %strcmp(ChannelL{i},'Lineage')
                c = [c i];
            end
        end
    end
    coef = [];
    sox2mask = AggrResults_primary.MeanNucSign(:,15) > 50;
    for i = 1:length(c)
        for j = 1:length(c)
            Ai = c(i);
            Bi = c(j);
            A = AggrResults_primary.MeanFullCellSign(:,Ai);
            B = AggrResults_primary.MeanFullCellSign(:,Bi);
            clearnegmask = A>0 & B>0 & sox2mask;% & soxmask;
            A = A(clearnegmask);
            B = B(clearnegmask);
            A = log10(A);
            B = log10(B);
            matrixab = [A,B];
            ccoef =corr(matrixab,'Type','Spearman');
            %ccoef = corrcoef(matrixab);
            coef{i+1,j+1} = round(ccoef(2),2,"significant");
        end
    end
    coef(2:2+length(c)-1,1) = AllChannel(c);
    coef(1,2:2+length(c)-1) = AllChannel(c);
    
    eucD = pdist(cell2mat(coef(2:end,2:end)),'euclidean');
    clustTreeEuc = linkage(eucD,'average');
    [h,nodes,outperm] = dendrogram(clustTreeEuc,0);
    
    c = c(outperm);
    coef = [];
    for i = 1:length(c)
        for j = 1:length(c)
            Ai = c(i);
            Bi = c(j);
            A = AggrResults_primary.MeanFullCellSign(:,Ai);
            B = AggrResults_primary.MeanFullCellSign(:,Bi);

            clearnegmask = A>0 & B>0 & sox2mask;% & soxmask &nonimmune;
            A = A(clearnegmask);
            B = B(clearnegmask);
            A = log10(A);
            B = log10(B);
            matrixab = [A,B];
            ccoef =corr(matrixab,'Type','Spearman');
            coef{i+1,j+1} = round(ccoef(2),2,"significant");
        end
    end
    coef(2:2+length(c)-1,1) = AllChannel(c);
    coef(1,2:2+length(c)-1) = AllChannel(c);
    
    %normalize
    cmat = cell2mat(coef(2:end,2:end));
    %Heatmap
    figure,
    h = heatmap(coef(1,2:end),coef(2:end,1)',cell2mat(coef(2:end,2:end)));
    set(gcf,'position',[400,100,1000,800])
    colormap('Hot')
    exportgraphics(gcf,[filename.graphsfolder 'Heatmap.jpg'],'Resolution',600)

end

