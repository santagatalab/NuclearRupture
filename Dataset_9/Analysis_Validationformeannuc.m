% Processing Multiple slides
clear all
filename.tissues = {'LSP15550', 'LSP15558'...
                    };    

for ctissuename_index = 1:length(filename.tissues)
    % 1) THE FOLDER DIRECTORIES
    ctissuename = filename.tissues{ctissuename_index}; % Current Tissue Name
    filename.tissue = ctissuename; 
    filename.datafolder = ['Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Liposarcoma\' ctissuename '\']; %Main Analysis location
    filename.rawfolder =  ['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Liposarcoma\' ctissuename '\']; %Ometif raw file location
    filename.graphsfolder = 'Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_02_Liposarcoma\Graphs\'; % Output Graph Folder
    filename.suffix = '.tif';
    filename.primaryNucfolder = [filename.datafolder 'ANALYSIS20230830\']; % Primary Nucleus Analysis folder
    filename.microNucfolder = [filename.datafolder 'ANALYSIS20230830micro\']; % Micronucleus Analysis folder
    filename.resultfolder = [filename.datafolder 'Results\']; % Result folder
    filename.fullstackfolder = [filename.datafolder 'FullStacks' filesep]; % Full Stacks Folder
    filename.folders.fullstacks = 'FullStacks\';
    filename.folders.coordinates = 'Coordinates_FullStacks\';
    filename.dim = 2; % Data Dimension
    
    % 2) SLIDE SPECIFIC PARAMETERS:
    if exist(filename.fullstackfolder,"dir")
        tempdir = dir(filename.fullstackfolder);
        [~,ind]=sort([tempdir.datenum]);
        tempdir = tempdir(ind);
        coren = 0;
        for i1 = 1:(length(dir(filename.fullstackfolder)))
            if contains(tempdir(i1).name,filename.tissues{ctissuename_index})
                coren = coren + 1;
                filename.realcoreinfo(coren).name = tempdir(i1).name;
                filename.realcoreinfo(coren).index = coren;
            end
        end
        clear i1 tempdir coren ind
    end
    filename.cycles = 12;  % cycles to analyse
    filename.maxround = filename.cycles*4;
    filename.ilastiksegfol = 'IlastikSegmentation\'; % Ilastik Folder setup
    filename.ilastiksuffix = '_Probabilities.tif';
    filename.segsuffix = '_Seg.tif';                    
    filename.ilastikround = 12; % Rounds to be used to segment in Ilastik (<= maxround)
    
    % 3) USER DESIRED PARAMETERS 
    filename.sizefield = 5000; %Size of field desired for a square 
    filename.crops = 1; %# of cropped fields desired per field 
    options.crops.DAPIbkgd_crops = 100;
    options.crops.prctile = 75;
    options.DAPIbkgd = 100; % thresh for DAPI used to check if we keep the field
    options.DAPIbkgd_crops75prctile = options.DAPIbkgd*5;
    
    % 4) OPTIONS  
    % Step 3: SEGMENTATION OPTIONS 
    options.nuc = 1;
    options.cyt = 2;
    options.backgr = 3;
    options.max_prob = 65535;
    options.bkgdprob_min = 0.2;
    options.cytprob_min  = 0.9;
    options.nucprob_min  = 0.7;
    options.nucprob_min_micro  = 0.6;
    options.pausetime = 0;
    
    
    % Step 4: MEASUREMENTS AND FOCI OPTIONS 
    options.tilesize = 1024;
    options.scale = 8;
    options.buffer = 256;
    options.coreradius = 200;
    options.dim = ['%02d'];
    options.rows = 2;
    options.cols = 2;
    options.background = 500;
    options.sigma = 0.3;
    options.alpha = 0.1;
    options.fileordchange = [1 2 3 4];
    options.cellsize = 100;
    options.cytosize = options.cellsize/5;
    options.corrthresh = 0.5;
    options.appendflag = 0;
    options.dates = '20230830_';
    
    filename.resultfile = [options.dates 'Results.mat'];
    filename.resultfilemicro = [options.dates 'microResults.mat'];
    
    
    options.substructures.flag = 1;
    options.substructures.channel = 20;
    options.substructures.bkgsize = 15;
    options.substructures.SV_thresh = 1;
    options.substructures.nucdilation = 3; % in case the nuclear segmentation is very stringent
    
    % 5) Cylinter
    filename.cylinterfol = 'Cylinter\';
    filename.cylinterinputfol = 'Input\';
    filename.cylintermaskfol = 'mask';
    filename.cylintersegfol = 'seg';
    filename.cylintertiffol = 'tif';
    
    % 6) Roi from Omero
    filename.roifile = [filename.datafolder 'Batch_ROI_Export'];
    
    % 7) Aggregate DATA
    options.cgaschannel = 42;

    load([filename.resultfolder 'Aggregate_' filename.resultfile])
    temp{ctissuename_index} = AggrResults_primary; 
    temp_MN{ctissuename_index} = AggrResults_micro; 
end

Channelcutoffcell = readcell([filename.datafolder 'Channelscutoff_BC.xlsx']);
[maxchannel,~] = size(Channelcutoffcell);
AllChannel = Channelcutoffcell(:,1);
AllChanneln = 1:maxchannel;
Channelcutoff = cell2mat(Channelcutoffcell(:,2));
amps = cell2mat(Channelcutoffcell(:,3));
%Channelchar = Channelcutoffcell(:,4);
ChannelL = Channelcutoffcell(:,5);
%% Thresholding
Core = 2;
cn = 6;

% Output an image. Red is positive foci, yellow is foci that is above
% threshhold but eliminated due to high surrounding signal.
Channelname =  AllChannel(cn);
Channel = AllChanneln(cn);
Threshold = 300;%Channelcutoff(cn);
amp = amps(cn);
mask = (AggrResults_primary.Corenum == Core) & (AggrResults_primary.MeanNucSign(:,Channel) > Threshold);
CentroidXlist = AggrResults_primary.CentroidX(mask) ;
CentroidYlist = AggrResults_primary.CentroidY(mask) ;

mask = (AggrResults_primary.Corenum == Core) & (AggrResults_primary.MeanNucSign(:,Channel) < Threshold);
CentroidXlist1 = AggrResults_primary.CentroidX(mask) ;
CentroidYlist1 = AggrResults_primary.CentroidY(mask) ;
% mask = (AggrResults.Name(:) == Core) & (AggrResults.MeanNucSign(:,Channel) > Threshold) & (AggrResults.AreaSubstruct(:,1)./AggrResults.AreaSubstruct(:,2) > 0.001);
% CentroidXlist = AggrResults.CentroidX(mask) ;
% CentroidYlist = AggrResults.CentroidY(mask) ;
% 
% mask = (AggrResults.Name(:) == Core) & (AggrResults.MeanNucSign(:,Channel) < Threshold) & (AggrResults.AreaSubstruct(:,1)./AggrResults.AreaSubstruct(:,2) > 0.001);
% CentroidXlist1 = AggrResults.CentroidX(mask) ;
% CentroidYlist1 = AggrResults.CentroidY(mask) ;


% %Threshold graphing
% edges = linspace(0, 500, 500); % Create 20 bins.
% % Plot the histogram.
% data = AggrResults.MeanFociSign(:,Channel)-double(AggrResults.MedianFociSign(:,Channel));
% data(data==0) = [];
% histogram(data, 'BinEdges',edges);
% % Fancy up the graph.
% grid on;
% xlim([0,100]);
% ylim([0,20000])
% xlabel('Data Value', 'FontSize', 14);
% ylabel('Bin Count', 'FontSize', 14);
% title('CGAS', 'FontSize', 14);
k1 = Core; %Core number 
Channel = Channel; %Interest Channel 
i3 = 12; %BANF Channel



%Start
fullstackfolder = filename.fullstackfolder;
FileTif = [fullstackfolder filename.realcoreinfo(k1).name];
DAPI_img = uint16(imread(FileTif,'Index',1));
ilastikfolder = [filename.primaryNucfolder filename.ilastiksegfol];
ilastiksegfolder = [ilastikfolder 'seg\'];
ilastikNucCytSegfolder = [ilastikfolder 'nuccytseg\'];
segfile = [ilastiksegfolder num2str(filename.realcoreinfo(k1).index) filename.segsuffix];
cyt_segfile = [ilastikNucCytSegfolder num2str(filename.realcoreinfo(k1).index) '_NucCytSeg.tif'];
NucMask = uint16(imread(segfile));
lb_NucMask = uint16(bwlabel(NucMask));

% create cytoplasmic mask
lb_CytMask = imdilate(lb_NucMask,offsetstrel('ball',ceil(options.cytosize),0));
lb_CytMask(lb_NucMask>0)=0;
CytMask = uint16(lb_CytMask > 0);


DAPI_img = uint16(imread(FileTif,'Index',1));
NucCytSeg = cat(3,NucMask,DAPI_img,CytMask);
imwrite(NucCytSeg,cyt_segfile);

%%%% measure cell properties and assign them to results -------
stats_nuc = regionprops(lb_NucMask,'Area','Centroid','Solidity');
stats_cyt = regionprops(lb_CytMask,'Area');

Area = {stats_nuc.Area};
Solidity = {stats_nuc.Solidity};
Centroid = {stats_nuc.Centroid};
CentroidVect = cell2mat(Centroid);
CentroidMat = reshape(CentroidVect,2,length(CentroidVect)/2);
CytArea = {stats_cyt.Area};

totcells_field = length(Area);
currentcells = 0; %length(Results{k1}.Area);
startcell = currentcells + 1;
endcell = currentcells + totcells_field;

FluorImage = imread(FileTif,'Index',Channel);
FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
close all
% imshow(DAPI_img,[])
% figure
% imshow(FociSegImage,[])
DAPI_img1 = DAPI_img*0;
image = cat(3,DAPI_img1,DAPI_img1,DAPI_img1);
% mask = cat(3,NucMask<0,NucMask<0,NucMask>0);
% image(mask) = 60000;
%SigSegImage = FluorImage_BK > 600;
% mask = cat(3,SigSegImage<0,SigSegImage,SigSegImage<0);
% image(mask) = 60000;
image(:,:,2) = image(:,:,2) + FluorImage*20;


centers = [CentroidXlist1,CentroidYlist1];
radii = zeros(length(CentroidXlist1),1) + 30;

% imshow(image)
%h = viscircles(centers,radii,'Color','r','LineWidth',.5);
radii = zeros(length(CentroidXlist1),1) + 30;
for i=1:length(radii)
    theta=0:1:360;
    r=round(centers(i,1) + radii(i)*sin(theta));
    c=round(centers(i,2) + radii(i)*cos(theta));
      for j=1:length(r)
          try
            image(c(j),r(j),:)=[60000 60000 0];
          end
      end
end
radii = zeros(length(CentroidXlist1),1) + 31;
for i=1:length(radii)
    theta=0:1:360;
    r=round(centers(i,1) + radii(i)*sin(theta));
    c=round(centers(i,2) + radii(i)*cos(theta));
      for j=1:length(r)
          try
            image(c(j),r(j),:)=[60000 60000 0];
          end
      end
end
centers = [CentroidXlist,CentroidYlist];
radii = zeros(length(CentroidXlist),1) + 30;

% imshow(image)
%h = viscircles(centers,radii,'Color','r','LineWidth',.5);
radii = zeros(length(CentroidXlist),1) + 30;
for i=1:length(radii)
    theta=0:1:360;
    r=round(centers(i,1) + radii(i)*sin(theta));
    c=round(centers(i,2) + radii(i)*cos(theta));
      for j=1:length(r)
          try
            image(c(j),r(j),:)=[60000 0 0];
          end
      end
end
radii = zeros(length(CentroidXlist),1) + 31;
for i=1:length(radii)
    theta=0:1:360;
    r=round(centers(i,1) + radii(i)*sin(theta));
    c=round(centers(i,2) + radii(i)*cos(theta));
      for j=1:length(r)
          try
            image(c(j),r(j),:)=[60000 0 0];
          end
      end
end
 figure,imshow(image)

path = ['Z:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\SavedimagepH2AX\' num2str(Core) '.tiff'];
options.append = true;
%saveastiff(image, path, options)
%end
%% Thresholding
Core = 57;
cn =20;


% Output an image. Red is positive foci, yellow is foci that is above
% threshhold but eliminated due to high surrounding signal.
Channelname =  AllChannel(cn);
Channel = AllChanneln(cn);
Threshold = 100;%Channelcutoff(cn);
amp = amps(cn);
mask = (AggrResults_micro.Name == Core) & (AggrResults_micro.MeanNucSign(:,Channel) > Threshold)  & (AggrResults_micro.MeanNucSign(:,Channel) > 1.3*AggrResults_micro.MeanCytSign(:,Channel));
CentroidXlist = AggrResults_micro.CentroidX(mask) ;
CentroidYlist = AggrResults_micro.CentroidY(mask) ;

mask = (AggrResults_micro.Name == Core) & (AggrResults_micro.MeanNucSign(:,Channel) > 0) & (AggrResults_micro.MeanNucSign(:,Channel) > 0*AggrResults_micro.MeanCytSign(:,Channel));
CentroidXlist1 = AggrResults_micro.CentroidX(mask) ;
CentroidYlist1 = AggrResults_micro.CentroidY(mask) ;
% mask = (AggrResults.Name(:) == Core) & (AggrResults.MeanNucSign(:,Channel) > Threshold) & (AggrResults.AreaSubstruct(:,1)./AggrResults.AreaSubstruct(:,2) > 0.001);
% CentroidXlist = AggrResults.CentroidX(mask) ;
% CentroidYlist = AggrResults.CentroidY(mask) ;
% 
% mask = (AggrResults.Name(:) == Core) & (AggrResults.MeanNucSign(:,Channel) < Threshold) & (AggrResults.AreaSubstruct(:,1)./AggrResults.AreaSubstruct(:,2) > 0.001);
% CentroidXlist1 = AggrResults.CentroidX(mask) ;
% CentroidYlist1 = AggrResults.CentroidY(mask) ;


% %Threshold graphing
% edges = linspace(0, 500, 500); % Create 20 bins.
% % Plot the histogram.
% data = AggrResults.MeanFociSign(:,Channel)-double(AggrResults.MedianFociSign(:,Channel));
% data(data==0) = [];
% histogram(data, 'BinEdges',edges);
% % Fancy up the graph.
% grid on;
% xlim([0,100]);
% ylim([0,20000])
% xlabel('Data Value', 'FontSize', 14);
% ylabel('Bin Count', 'FontSize', 14);
% title('CGAS', 'FontSize', 14);
k1 = Core; %Core number 
Channel = Channel; %Interest Channel 
i3 = 12; %BANF Channel



%Start
fullstackfolder = [filename.datafolder 'dearray' filesep]; 
FileTif = [fullstackfolder num2str(filename.realcoreinfo(k1).index)  '.ome.tif'];
DAPI_img = uint16(imread(FileTif,'Index',1));
ilastikfolder = [filename.analfolder filename.ilastiksegfol];
ilastiksegfolder = [ilastikfolder 'seg\'];
ilastikNucCytSegfolder = [ilastikfolder 'nuccytseg\'];
segfile = [ilastiksegfolder num2str(filename.realcoreinfo(k1).index) filename.segsuffix];
cyt_segfile = [ilastikNucCytSegfolder num2str(filename.realcoreinfo(k1).index) '_NucCytSeg.tif'];
NucMask = uint16(imread(segfile));
lb_NucMask = uint16(bwlabel(NucMask));

% create cytoplasmic mask
lb_CytMask = imdilate(lb_NucMask,offsetstrel('ball',ceil(options.cytosize),0));
lb_CytMask(lb_NucMask>0)=0;
CytMask = uint16(lb_CytMask > 0);


DAPI_img = uint16(imread(FileTif,'Index',1));
NucCytSeg = cat(3,NucMask,DAPI_img,CytMask);
imwrite(NucCytSeg,cyt_segfile);

%%%% measure cell properties and assign them to results -------
stats_nuc = regionprops(lb_NucMask,'Area','Centroid','Solidity');
stats_cyt = regionprops(lb_CytMask,'Area');

Area = {stats_nuc.Area};
Solidity = {stats_nuc.Solidity};
Centroid = {stats_nuc.Centroid};
CentroidVect = cell2mat(Centroid);
CentroidMat = reshape(CentroidVect,2,length(CentroidVect)/2);
CytArea = {stats_cyt.Area};

totcells_field = length(Area);
currentcells = 0; %length(Results{k1}.Area);
startcell = currentcells + 1;
endcell = currentcells + totcells_field;

FluorImage = imread(FileTif,'Index',Channel);
FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
% imshow(DAPI_img,[])
% figure
% imshow(FociSegImage,[])
DAPI_img1 = DAPI_img*3;
image = cat(3,DAPI_img1,DAPI_img1,DAPI_img1);
mask = cat(3,NucMask<0,NucMask<0,NucMask>0);
%image(mask) = 60000;
SigSegImage = FluorImage_BK > 600;
% mask = cat(3,SigSegImage<0,SigSegImage,SigSegImage<0);
% image(mask) = 60000;
% image(:,:,1) = image(:,:,2) + FluorImage*50;

%BAF
FluorImage = imread(FileTif,'Index',63);
FluorImage_BK = FluorImage - imopen(imclose(FluorImage,strel('disk',ceil(options.cellsize/10))),strel('disk',ceil(options.cellsize*3)));
image(:,:,2) = image(:,:,2) + FluorImage*3;


centers = [CentroidXlist1,CentroidYlist1];
radii = zeros(length(CentroidXlist1),1) + 30;

% imshow(image)
%h = viscircles(centers,radii,'Color','r','LineWidth',.5);
radii = zeros(length(CentroidXlist1),1) + 30;
for i=1:length(radii)
    theta=0:1:360;
    r=round(centers(i,1) + radii(i)*sin(theta));
    c=round(centers(i,2) + radii(i)*cos(theta));
      for j=1:length(r)
          try
            image(c(j),r(j),:)=[60000 60000 0];
          end
      end
end
radii = zeros(length(CentroidXlist1),1) + 31;
for i=1:length(radii)
    theta=0:1:360;
    r=round(centers(i,1) + radii(i)*sin(theta));
    c=round(centers(i,2) + radii(i)*cos(theta));
      for j=1:length(r)
          try
            image(c(j),r(j),:)=[60000 60000 0];
          end
      end
end
centers = [CentroidXlist,CentroidYlist];
radii = zeros(length(CentroidXlist),1) + 30;

% imshow(image)
%h = viscircles(centers,radii,'Color','r','LineWidth',.5);
radii = zeros(length(CentroidXlist),1) + 30;
for i=1:length(radii)
    theta=0:1:360;
    r=round(centers(i,1) + radii(i)*sin(theta));
    c=round(centers(i,2) + radii(i)*cos(theta));
      for j=1:length(r)
          try
            image(c(j),r(j),:)=[60000 0 0];
          end
      end
end
radii = zeros(length(CentroidXlist),1) + 31;
for i=1:length(radii)
    theta=0:1:360;
    r=round(centers(i,1) + radii(i)*sin(theta));
    c=round(centers(i,2) + radii(i)*cos(theta));
      for j=1:length(r)
          try
            image(c(j),r(j),:)=[60000 0 0];
          end
      end
end
 figure,imshow(image)

path = ['Z:\sorger\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\ANALYSIS May 2022\SavedimagepH2AX\' num2str(Core) '.tiff'];
options.append = true;
%saveastiff(image, path, options)
%end
%%
c = [];
for i = 1:length(Channelchar)
    try
        if (strcmp(Channelchar{i},'Nuclear') || strcmp(Channelchar{i},'Cytoplasmic')) && ChannelL{i} == 'L' %Palentirmarks{i} == 'Palantir' %
            c = [c i];
        end
    end
end
for i = 17:23
    subplot(8,1,i-16)
    hold on
    ccn = AllChanneln(c(i));
    aggrcurrent = AggrResults_primary.MeanFullCellSign(:,ccn);
    ksdensity(log10(aggrcurrent+1))
    xlabel(['log10(Mean Nuclear Magnitude of ' AllChannel{c(i)} ')'])
    if i == 1
       title('Probability Density Distribution Comparison between Micronucleus and Primary Nucleus')
    end
    hold off
end

