% stage3
%
%   stage3 runs the third stage of the kilosort spike sorting for an
%   interactive job on Big Red 2. This stage handles imposing the rulings 
%   from stage 2 on the data, organizing the data, and shipping it back to 
%   IU Box.

%% Settings

% Set the main directory on the data capacitor to store data (something
% like...)
mainDC = '/N/dc2/scratch/nmtimme/ver2IntTest1';

%% Perform Necessary Tasks

% Load the stage 1 job information
load([mainDC,filesep,'spikeSortingStage1Info.mat'])

%% Make a file with information 

% Save the file in your main Karst directory so it is easy to find
save([mainDC,filesep,'spikeSortingStage3Info.mat'],'boxDataSetDir','dataSetParams','dcDataSetDir','dataSetID','IUstring')

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
stage3Core(info)