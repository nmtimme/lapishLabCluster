% stage1JobsToDo
%
%   stage1JobsToDo creates a list of KiloSort spike sorting 
%   jobs to perform on a computer cluster for stage 1 of the analysis. 
%   This stage handles moving data from IU Box to the Data Capacitor and 
%   generating a report on basic features of the data to be reviewed by 
%   the user.

%% Settings

% Set the main directory on the data capacitor to store data (something
% like...)
mainDC = '/N/dc2/scratch/nmtimme/ver2ParTest1';


% Set the queue that should receive the jobs (leave empty if not requesting
% a certain queue)
queueName = [];
% queueName = '-q debug '; % debug queue
% queueName = '-q preempt '; % preempt queue (only on Karst)
% queueName = '-q gpu '; % GPU enabled nodes (only on Big Red 2)


%% Data Set Information

% Note: the values in boxDataSetDirs, dataSetIDs, and dataSetParams are
% linked such that boxDataSetDirs{i}, dataSetIDs{i}, and dataSetParams(i,:)
% all correspond to the same spike sorting analysis. To run spike sorting
% on multiple analyses, increase the length of the these arrays with each
% new spike sorting analysis.

% List the directories of each data set on box (identical data sets could
% be listed more than once with different parameters, see below)
boxDataSetDirs = {'NickEphysData/2018-11-12 - Session3/openEphysData';...
    'NickEphysData/2018-11-12 - Session3/openEphysData'};



% List the names of the data sets (identical data sets names could be 
% listed more than once with different parameters, see below) (note, no
% spaces should be used)
dataSetIDs = {'2018-11-12-Session3ParNoRef';...
    '2018-11-12-Session3ParAveRef'};



% Set parameters for each data set. This will be stored as an N by M cell
% array (N = number of data sets, M = number of parameters). Note, many 
% other parameters are hard coded into the stage 2 job script. Here is the
% meaning of each parameter:
% dataSetParams(:,1) = number of channels in the whole recording (usually 16, 32, or 64)
% dataSetParams(:,2) = parameter that influences the number of clusters per channel for initial sort (roughly the number of clusters per channel)
% dataSetParams(:,3) = spike threshold in standard deviations (we aren't sure if this is even used)
% dataSetParams(:,4) = referencing option (0: no referencing, 1: common average referencing, 2: common median referencing (too slow to use))
% dataSetParams(:,5) = maximum number of spikes to extract per batch (20000)
% dataSetParams(:,6) = threshold for detecting spikes on template-filtered data ([4 10 10])
% dataSetParams(:,7) = large means amplitudes are forced around the mean ([5 20 20])
% 
dataSetParams = {64,4,-3,0,20000,[4 10 10],[5 20 20];...
    64,4,-3,1,20000,[4 10 10],[5 20 20]};


% Alert the user to the number of jobs
nJobs = length(boxDataSetDirs);
disp(['This batch consists of ',num2str(nJobs),' jobs.'])

%% Make directories for the data

% Throw an error if the main data capacitor directory already exists
if exist(mainDC,'dir') == 7
    error('Data Capacitor directory already exists. Delete or rename.')
else
    mkdir(mainDC)
end

% Make the directories on the data capacitor for each data set
dcDataSetDirs = cell(size(boxDataSetDirs));
for iJob = 1:nJobs
    dataSetDir = [mainDC,filesep,num2str(iJob),'-',dataSetIDs{iJob}];
    mkdir(dataSetDir)
    dcDataSetDirs{iJob} = dataSetDir;
end

% Make directories on the data capacitor to save the pre and post user
% reviews of the results of the stage 1 jobs
mkdir([mainDC,'/Stage1ResultsPreReview'])
mkdir([mainDC,'/Stage1ResultsPostReview'])

%% Get IU Box Password

% Request that the user enter the password for IU Box
IUstring = input('Enter IU Box Password: ','s');


%% Make a file with information for the jobs

% Save the file in your main Karst directory so it is easy to find
save([mainDC,filesep,'spikeSortingStage1Info.mat'],'boxDataSetDirs','dataSetParams','dcDataSetDirs','dataSetIDs','IUstring')



%% Make a text file to start all the jobs on Karst

% Go back to the main director
matDir = pwd;
cd ~

% Make the text cell
TxtCell = cell(nJobs + 1,1);

% Put in the header
TxtCell{1,1} = '#!/bin/bash';

% List all the jobs
for iJob = 1:nJobs
    
    TxtCell{iJob + 1,1} = ['qsub ',queueName,'-t ',num2str(iJob),' spikeSortStage1JobVer1.txt']; % NOTE: THE FILE NAME MUST MATCH YOUR JOB TEXT FILE
    
end

% Name the job list
txtfilename = 'spikeSortStage1JobList.txt'; % NOTE: THIS WILL SET THE NAME OF YOUR JOB LIST (SEE NOTE BELOW TO EXECUTE THE LIST)

% Write the information to the text file
fid = fopen(txtfilename, 'w');
for iJob = 1:(nJobs + 1)
    fprintf(fid, '%s\r\n', TxtCell{iJob,:});
end
fprintf(fid,'\r\n');
fclose(fid);


% Note, to start all the jobs, enter the following commands in your main
% Big Red II directory
% chmod u+x spikeSortStage1JobList.txt
% dos2unix spikeSortStage1JobList.txt
% ./spikeSortStage1JobList.txt

cd(matDir)
