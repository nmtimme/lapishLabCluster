% stage1Job
%
%   stage1Job runs the first stage of the kilosort spike sorting on Big 
%   Red 2 for multiple parallel data sets. This stage handles moving data 
%   from IU Box to the Data Capacitor and generating a report on basic 
%   features of the data to be reviewed by the user.

%% Settings

% List any general settings for the spike sorting here. These won't change
% between jobs.

% Set the main directory on the data capacitor to store data (something
% like...)
mainDC = '/N/dc2/scratch/edefalco/batch6';


%% Load the job information

% Find the job number
jobNum = str2num(getenv('PBS_ARRAYID'));
% jobNum = 1;

% Load the job info
load([mainDC,filesep,'spikeSortingStage1Info.mat'])

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
info.username = username;
info.suffix= suffix;

% Run the core function
stage1Core(info)


