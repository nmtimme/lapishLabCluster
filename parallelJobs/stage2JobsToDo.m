% stage2JobsToDo
%
%   stage2JobsToDo creates a list of KiloSort spike sorting jobs to 
%   perform on a computer cluster for stage 2 of the analysis. This stage 
%   handles preparing data for spike sorting, spike sorting, and
%   generating a report of the spike sorting results.


%% Settings

% Set the main directory on the data capacitor to store data (something
% like...)
mainDC = '/N/dc2/scratch/edefalco/batch6';

% Set the queue that should receive the jobs (leave empty if not requesting
% a certain queue)
% queueName = [];
% queueName = '-q debug '; % debug queue
% queueName = '-q preempt '; % preempt queue (only on Karst)
queueName = '-q gpu '; % GPU enabled nodes (only on Big Red 2)

%% Perform Necessary Tasks

% Load the stage 1 job information
load([mainDC,filesep,'spikeSortingStage1Info.mat'])

% Alert the user to the number of jobs
nJobs = length(boxDataSetDirs);
disp(['This batch consists of ',num2str(nJobs),' jobs.'])

% Make directories on the data capacitor to save the pre and post user
% reviews of the results of the stage 2 jobs
mkdir([mainDC,'/Stage2ResultsPreReview'])
mkdir([mainDC,'/Stage2ResultsPostReview'])

%% Make a file with information for the jobs

% Save the file in your main Karst directory so it is easy to find
save([mainDC,filesep,'spikeSortingStage2Info.mat'],'boxDataSetDirs','dataSetParams','dcDataSetDirs','dataSetIDs','IUstring','username','suffixs')

%% Make a text file to start all the jobs on Karst

% Go back to the main directory
matDir = pwd;
cd ~

% Make the text cell
TxtCell = cell(nJobs + 1,1);

% Put in the header
TxtCell{1,1} = '#!/bin/bash';

% List all the jobs
for iJob = 1:nJobs
    
    TxtCell{iJob + 1,1} = ['qsub ',queueName,'-t ',num2str(iJob),' spikeSortStage2JobVer1.txt']; % NOTE: THE FILE NAME MUST MATCH YOUR JOB TEXT FILE
    
end

% Name the job list
txtfilename = 'spikeSortStage2JobList.txt'; % NOTE: THIS WILL SET THE NAME OF YOUR JOB LIST (SEE NOTE BELOW TO EXECUTE THE LIST)

% Write the information to the text file
fid = fopen(txtfilename, 'w');
for iJob = 1:(nJobs + 1)
    fprintf(fid, '%s\r\n', TxtCell{iJob,:});
end
fprintf(fid,'\r\n');
fclose(fid);


% Note, to start all the jobs, enter the following commands in your main
% Karst directory
% chmod u+x spikeSortStage2JobList.txt
% dos2unix spikeSortStage2JobList.txt
% ./spikeSortStage2JobList.txt

% Go back to the matlab directory
cd(matDir)
