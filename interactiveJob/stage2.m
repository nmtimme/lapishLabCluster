% stage2
%
%   stage2 runs the second stage of the kilosort spike sorting for an
%   interactive job on Big Red 2. This stage handles preparing data for 
%   spike sorting, spike sorting, and generating a report of the spike 
%   sorting results.

%% Settings

% Set the main directory on the data capacitor to store data (something
% like...)
mainDC = '/N/dc2/scratch/nmtimme/ver2IntTest2';

%% Perform Necessary Tasks

% Load the stage 1 job information
load([mainDC,filesep,'spikeSortingStage1Info.mat'])

% Make directories on the data capacitor to save the pre and post user
% reviews of the results of the stage 2 jobs
mkdir([mainDC,'/Stage2ResultsPreReview'])
mkdir([mainDC,'/Stage2ResultsPostReview'])

%% Make a file with information 

% Save the file in your main Karst directory so it is easy to find
save([mainDC,filesep,'spikeSortingStage2Info.mat'],'boxDataSetDir','dataSetParams','dcDataSetDir','dataSetID','IUstring','username')

%% Run the core function

% Add the clustering software directories to the path
cd ..
addpath(genpath(pwd))

% Put the necessary information in an input structure for the core function
info = struct;
info.boxDataSetDir = boxDataSetDir;
info.dataSetParams = dataSetParams;
info.dcDataSetDir = dcDataSetDir;
info.dataSetID = dataSetID;
info.IUstring = IUstring;
info.mainDC = mainDC;

% Run the core function
stage2Core(info)