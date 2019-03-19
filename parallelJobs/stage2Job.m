% stage2Job
%
%   stage2Job runs the second stage of the kilosort spike sorting on Big 
%   Red 2 for multiple parallel data sets. This stage handles preparing 
%   data for spike sorting, spike sorting, and generating a report of the 
%   spike sorting results.

%% Settings

% List any general settings for the spike sorting here. These won't change
% between jobs.

% Set the main directory on the data capacitor to store data (something
% like...)
mainDC = '/N/dc2/scratch/edefalco/batch6';


%% Load the job information

% Find the job number
jobNum = str2num(getenv('PBS_ARRAYID'));
%jobNum = 16;

% Load the job info
load([mainDC,filesep,'spikeSortingStage2Info.mat'])

% Get the directories and parameters for this specific job
boxDataSetDir = boxDataSetDirs{jobNum};
params = dataSetParams(jobNum,:);
dcDataSetDir = dcDataSetDirs{jobNum};
dataSetID = dataSetIDs{jobNum};
suffix =suffixs{jobNum};


%% Run the core function

% Add the clustering software directories to the path
cd ..
addpath(genpath(pwd))

% Put the necessary information in an input structure for the core function
info = struct;
info.boxDataSetDir = boxDataSetDir;
info.dataSetParams = params;
info.dcDataSetDir = dcDataSetDir;
info.dataSetID = dataSetID;
info.IUstring = IUstring;
info.mainDC = mainDC;
info.suffix=suffix;

% Run the core function
stage2Core(info)

