% stage1
%
%   stage1 runs the first stage of the kilosort spike sorting for an
%   interactive job on Big Red 2. This stage handles moving data from IU
%   Box to the Data Capacitor and generating a report on basic features of
%   the data to be reviewed by the user.

%% Settings

% Set the main directory on the data capacitor to store data (something
% like...)
mainDC = '/N/dc2/scratch/nmtimme/ver2IntTest1';


%% Data Set Information

% List the directory of the data set on box 
boxDataSetDir = 'NickEphysData/2018-11-12 - Session3/openEphysData';

% List the name of the data set (note, no
% spaces should be used)
dataSetID = '2018-11-12-Session3';

% Set parameters for the data set. This will be stored as an 1 by M cell
% array (M = number of parameters). Note, many other parameters are hard 
% coded into the stage 2 job script. Here is the meaning of each parameter:
% dataSetParams(1) = number of channels in the whole recording (usually 16, 32, or 64)
% dataSetParams(2) = parameter that influences the number of clusters per channel for initial sort (roughly the number of clusters per channel)
% dataSetParams(3) = spike threshold in standard deviations (we aren't sure if this is even used)
% dataSetParams(4) = referencing option (0: no referencing, 1: common average referencing, 2: common median referencing)
% dataSetParams(5) = maximum number of spikes to extract per batch (20000)
% dataSetParams(6) = threshold for detecting spikes on template-filtered data ([4 10 10])
% dataSetParams(7) = large means amplitudes are forced around the mean ([5 20 20])
% 
dataSetParams = {64,4,-3,1,20000,[4 10 10],[5 20 20]};



%% Make directories for the data

% Throw an error if the main data capacitor directory already exists
if exist(mainDC,'dir') == 7
    error('Data Capacitor directory already exists. Delete or rename.')
else
    mkdir(mainDC)
end

% Make the directories on the data capacitor for each data set
dataSetDir = [mainDC,filesep,dataSetID];
mkdir(dataSetDir)
dcDataSetDir = dataSetDir;

% Make directories on the data capacitor to save the pre and post user
% reviews of the results of the stage 1 jobs
mkdir([mainDC,'/Stage1ResultsPreReview'])
mkdir([mainDC,'/Stage1ResultsPostReview'])

%% Get IU Box Password

% Request that the user enter the password for IU Box
IUstring = input('Enter IU Box Password: ','s');


%% Make a file with information 

% Save the file in the main Data Capacitory directory for this data set so
% it is easy to find
save([mainDC,filesep,'spikeSortingStage1Info.mat'],'boxDataSetDir','dataSetParams','dcDataSetDir','dataSetID','IUstring')

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
stage1Core(info)

