% t-CycIF ASHLAR ANALYSIS PIPELINE 
% Edited by Brian Cheng 07/21/2023
%
% Options:  Primary Nucleus Only/Micronucleus 
%           Whole Slide Imaging (WSI)/Tissue Microarrays (TMA)
%
% Step 1: (Step1_fields_preilastik.m) Cut the "ome.tif" into 
%         fields of a size specified by the user and cut crops out of each field 
%         for Ilastik training. Create full stacks of all cycles and channels. 
%         Omits fields that have no cells.   
% Step 2: (Step2_filtercrops_v2.m) Remove blank regions 
% Step 3: (Step3_segmentfromilastik.m) Segment based on segmentation
%         probabilities produce by Ilastik
% Step 4: (Step4_CycIF_measurements_micro.m) Makes measurements of
%         signal and foci segmentation 
% Step 5: (Step5_Cylinter_Result.m) [Optional] Output Result to Cylinter
%         Create a folder for Ilastik analysis
%   then: (Step5_Cylinter_Resultout.m) [Optional] Input Result from Cylinter
% Step 6: (Step6_ROI_Omero.m) [Optional] Input ROI from Omero
%         Take ROI batch export from Omero and remove cells from Results
%         that is not in the ROI
% Step 7: (Step7_Analysis_Prep.m) Aggregate data 
%         Create AggregateResults and other data structure for downstream
%         analysis

% DATA Analysis and Graph Generation
% Graph Output Folder: <BaseFolder>\Graphs\*.pdf
% Option 1: (Analysis_CrossCorrelation.m) Create a marker Heatmap 
% Option 2: (Analysis_BAFandMicronucleusComparison.m) Comparison graph
%           between different groups: BAF+/BAF- MN+/MN- MNBAF+/MNBAF-
% Option 3: (Analysis_TumorBoundary.m) Spacial analysis on Tumor Boundary
%           Required: (p_poly_dist.m) 

% In order to begin the analysis a set of parameters are needed to be
% defined by the user 
% 1) where the files are located and how the names are formatted
% 2) parameters to choose the size of the field and number of crops per
% field 

%%% OUTPUTS FILES AND LOCATIONS 
% Step 1:
% FullStacks: <MainFolder>\FullStacks\OmeTifName_Field_row_column.tif
% Coordinates: <MainFolder>\Coordinates_FullStacks\OmeTifName.mat 
% y= # of pixels in 1 column; x = # of pixels in 1 row ; t = size of field desired 
% Coordinates Matrix: Coordinates.Field(row,column) = [keep field (0= no, 1=yes), x1, x2, y1, y2]
% Step 3:
% Segmented Images: <MainFolder>\ANALYSIS<date>\IlastikSegmentation\seg\*.tif
% Check Segmented Images: <MainFolder>\ANALYSIS\IlastikSegmentation\check\*.tif
% Step 4: 
% Foci Segmented Images: <MainFolder>\ANALYSIS<date>\IlastikSegmentation\nuccytseg\*.tif
% Foci Check Segmented Images: <MainFolder>\ANALYSIS<date>\IlastikSegmentation\checknucsubseg\*.tif
% Measurments file: <MainFolder>\Results\<date>_Results.mat 

%%% INPUT FILE LOCATION AND FORMATTING
%
% The code will run on the Ilastik Probabilities file and FullStacks file 
% The expected file structure is that a master folder will contain all of
% the "ome.tif" data and within that folder, a folder called "ANALYSIS\"
% will contain all the Ilastik Probabilities and FullStacks data 
%
%%% IMPORTANT: Ilastik Output Filename Format for Export Options
% AllTheRawData\ANALYSIS\MyProject.ilp 
% AllTheRawData\ANALYSIS\IlastikSegmentation\pmap\*Probabilities.tif 

% ROI:
% Can be created through Cylinter (Step 5) and/or Omero (Step 6)
% Omero batch export Location: AllTheRawData\Batch_ROI_Export.csv

% Processing Multiple slides
clear all
filename.tissues = {'LSP15749', 'LSP15756', 'LSP15762', 'LSP15767' ...
                    };    

for ctissuename_index = 1:length(filename.tissues)
    % 1) THE FOLDER DIRECTORIES
    ctissuename = filename.tissues{ctissuename_index}; % Current Tissue Name
    filename.tissue = ctissuename; 
    filename.datafolder = ['Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_04_GBM_Lineagev2\' ctissuename '\']; %Main Analysis location
    filename.rawfolder =  ['Y:\lsp-data\cycif-production\132-Nuclear-atypia-and-BAF\2023_03_GBM_TMA_Lineagev2\' ctissuename '\']; %Ometif raw file location
    filename.graphsfolder = 'Y:\lsp-analysis\cycif-production\132-Nuclear-atypia-and-BAF\2023_04_GBM_Lineagev2\Graphs\'; % Output Graph Folder
    filename.suffix = '.tif';
    filename.primaryNucfolder = [filename.datafolder 'ANALYSIS20230419\']; % Primary Nucleus Analysis folder
    filename.microNucfolder = [filename.datafolder 'ANALYSIS20230419micro\']; % Micronucleus Analysis folder
    filename.resultfolder = [filename.datafolder 'Results\']; % Result folder
    filename.fullstackfolder = [filename.datafolder 'FullStacks' filesep]; % Full Stacks Folder
    filename.folders.fullstacks = 'FullStacks\';
    filename.folders.coordinates = 'Coordinates_FullStacks\';
    filename.Channelcutoff = 'Channelscutoff_BC.xlsx';
    filename.dim = 2; % Data Dimension
    
    % 2) SLIDE SPECIFIC PARAMETERS:
    if exist(filename.fullstackfolder,"dir")
        tempdir = dir(filename.fullstackfolder);
        [~,ind]=sort([tempdir.datenum]);
        tempdir = tempdir(ind);
        coren = 0;
        for i1 = 1:(length(dir(filename.fullstackfolder)))
            if contains(tempdir(i1).name,filename.tissues{1})
                coren = coren + 1;
                filename.realcoreinfo(coren).name = tempdir(i1).name;
                filename.realcoreinfo(coren).index = coren;
            end
        end
        clear i1 tempdir coren ind
    end
    filename.cycles = 17;  % cycles to analyse
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
    options.dates = '20230720_';
    
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

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Recommend Process the following functions one-by-one and quality check
    %each step
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Data Processing functions
%     % If WSI, crop data using Step1 and Step 2
%     Step1_fields_preilastik(filename,options) 
%     Step2_filtercrops_v2(filename,options) 
%     % After running Step2, export FullStacks into Ilastik to train 
%     % If TMA, directely import TMA core fullstack data from MCMICRO 
%     % pipeline to Ilastik
%     % After obtaining probability maps, import files into filename.ilastiksegfol\pmap\
%     Step3_segmentfromilastik(filename,options) 
%     Step4_CycIF_measurements(filename,options) 
%     % After running Step4, automatically export to Cylinter format through
%     % Step5_Cylinter_Result. Import Cylinter output through Step5_Cylinter_Resultout
%     Step5_Cylinter_Result(filename,options)
%     Step5_Cylinter_Resultout(filename,options)
%     % Optional micronucleus processing before step6
%     Step3_segmentfromilastik_micro(filename,options) 
%     Step4_CycIF_measurements_micro(filename,options) 
%     Step5_Cylinter_Result_micro(filename,options)
%     Step5_Cylinter_Resultout_micro(filename,options)
%     % After running Step5, import roifile from Omero server as filename.roifile
%     Step6_ROI_Omero(filename,options)
%     Step7_Analysis_Prep(filename,options)
% 
%     % Generating Results Graphs
%     Analysis_CrossCorrelation(filename,options)
%     Analysis_BAFandMicronucleusComparison(filename,options)
%     Analysis_TumorBoundary(filename,options)
    
%     % For reproducibility, use the following function for past data
%     % Add dir: \Codes\DatasetSpecific

%     % Dataset_3 TMAs:
%     % Import DATA from MCMICRO
%       Dataset_3_Step3_segmentfromilastik_TMA(filename,options)
%       Dataset_3_Step4_CycIF_measurements_TMA(filename,options)
%       Dataset_3_Step5_ResulttoCylinter(filename,options)
%       Dataset_3_Step6_CylinterResultout(filename,options)
%     % Put the following files within the filename.analfolder
%     % Two columns with first row being the result file core number,
%     % second row being the TMA original number 
%       filename.CoreAssign = 'corenumbercorrection.xlsx';
%     % Channel Information: Channelname, Channelcutoff, Channelcutoff
%     % multiplier, Nuclear/Cytoplasmic signal, ChannelType
%       filename.Channelcutoff = 'Channelscutoff_BC.xlsx';
%     % TMA file information
%       filename.TMAclasslist = 'HTMA399_Diagnosis_Map_DeID_v4.xls';
%       Dataset_3_Analysis_LineageComparison(filename,options)

%     % Dataset_2 TMAs:
%     % Import DATA from MCMICRO
%       Dataset_2_Step3_segmentfromilastik_TMA(filename,options)
%       Dataset_2_Step4_CycIF_measurements_TMA(filename,options)
%       Dataset_2_Step5_ResulttoCylinter(filename,options)
%       Dataset_2_Step6_CylinterResultout(filename,options)
%     % Need to Adjust Individual Variables for Dataset_2
%       Dataset_2_Analysis_Bafanalysisaggregate

%     % Dataset_1 TMAs:
%     % Require Different Run-ALL
%     % Dataset_1_RUN_TMA_imageanalysis.m
%       Dataset_1_Step1_RealCoreCreation(filename,options)
%       Dataset_1_Step2_makecorestacks(filename,options)
%       Dataset_1_Step3_segmentfromilastik_TMA(filename,options)
%       Dataset_1_Step4_CycIF_measurements_TMA(filename,options)
%       Dataset_1_Step5_Analysis_Prep(filename,options)
%     % Need to Adjust Individual Variables for Dataset_1
%       Dataset_1_Analysis_GraphandTable

end