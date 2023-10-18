%This program runs the image analysis for TMAs
%Step1: Making the core montage using 1 10x simage- finds each core and
%outputs a RealCore tif file for each core
%Step 2: Makes the CoreStacks for each core: ouputs DAPI_Core (DAPIStack) &
%Core (FullStack)
%Step 3: Segments each Core: outputs TrackedField stack of segmented DAPI
%images 

% parameters that should not change
clear all


%Change All parameters below  


filename.suffix = '.tif';
filename.datafolder = 'Z:\data\IN_Cell_Analyzer_6000\Giorgio\2019-10-18 BANF1 GBM TMA Analysis\GBM_TMA_399_1_LONG_03_28_18_40X\Analysis\';
filename.analfolder = [filename.datafolder 'ANALYSIS May 2022\'];
filename.resultfolder = [filename.analfolder 'Results\'];
filename.codesfolder = [filename.analfolder 'Codes Used\'];

filename.wavelengths = {' wv UV - DAPI)',' wv Blue - FITC)',' wv Green - dsRed)',' wv Red - Cy5)'};
filename.prefix1 = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O'}; %
filename.prefix2 = linspace(1,25,25);
filename.dim = ['%02d']; 
filename.dim2 = ['%01d']; % this is for the number of tiles in a core
filename.cycleprefix = 'Cycle';
filename.realcorecycleprefix = 'Cycle1';
filename.cycles = [1 2 3 4 5 6 7];  % cycles to analyse
filename.maxcycle = 8;


filename.ilastiksegfol = 'Ilastik Segmentation\';
filename.ilastiksuffix = '_Probabilities.tif';
filename.segsuffix = '_Seg.tif';

options.nuc = 1;
options.cyt = 2;
options.backgr = 3;
options.max_prob = 65535;
options.bkgdprob_min = 0.2;
options.cytprob_min  = 0.8;
options.nucprob_min  = 0.33;
options.pausetime = 0;


count = 0;
for i1 = 1:length(filename.prefix1)
    for i2 = 1:length(filename.prefix2)
        count = count + 1;
        filename.realcoreinfo{count}.name = [filename.prefix1{i1} ' - ' num2str(filename.prefix2(i2),filename.dim) ];
        filename.realcoreinfo{count}.tiles = [1 2 3 4];
    end
end

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
options.cellsize = 35;
options.cytosize = options.cellsize/5;
options.corrthresh = 0.5;
options.appendflag = 0;
options.dates = '2023-1-20_';

filename.resultfile = [options.dates 'Results.mat'];


options.substructures.flag = 1;
options.substructures.channel = 12;
options.substructures.bkgsize = 15;
options.substructures.SV_thresh = 3;
options.substructures.nucdilation = 3; % in case the nuclear segmentation is very stringent

DAPIslice = (1:filename.maxcycle)*4 -3;

%For 10X Montage
% rows = 1;
% cols = 1;
% maxtile = rows*cols;
cut = 100;
background = 1000;
radius10X = 700;

% images per each core
rows2 = 2;
cols2 = 2;
maxtile2 = rows2*cols2;

% if this is only an update modify this parameters
num_old_DAPI = 0;


%% Running image Analysis

%Step1_RealCoreCreation(filename,options)
%%
%Step2_makecorestacks(filename,options)
%%
%Step3_segmentfromilastik_TMA(filename, options) 
% Step3_segmentfromilastik_TMA_label(filename, options) 
%%
%Step4_CycIF_measurements_TMA(filename, options) 
Step4_CycIF_measurements_TMA_label(filename, options) 


%% Testing 6/3
