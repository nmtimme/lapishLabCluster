function [] = stage3Core(info)
%STAGE3CORE runs the third stage of the kilosort spike sorting on Big Red 
%   2. This stage handles imposing the rulings from stage 2 on the data, 
%   organizing the data, and shipping it back to IU Box.
%
%   info: structure containing information about the data set to be spike
%   sorted.
%       info.boxDataSetDir (string): Full file path for the directory on IU
%           Box where the data set is stored.
%       info.dataSetParams (1xN cell array): Parameters for spike sorting.
%           See stage1JobsToDo.m or stage1.m.
%       info.dcDataSetDir (string): Full file path for the directory on the
%           Data Capacitor where the data will be stored.
%       info.dataSetID (string): Temporary name for the data set.
%       info.IUstring (string): IU Box Password.
%       info.mainDC (string): Full file path for the main spike sorting
%           directory on the Data Capacitor.

%% Pull variables out of the info structure

boxDataSetDir = info.boxDataSetDir;
dataSetParams = info.dataSetParams;
dcDataSetDir = info.dcDataSetDir;
dataSetID = info.dataSetID;
IUstring = info.IUstring;
mainDC = info.mainDC;

%% Load the results of the user review of the stage 2 results for this job

% Copy the post user review stage 2 results file into the directory for
% this specific analysis to ease transfer back to IU Box
copyfile([mainDC,filesep,'Stage2ResultsPostReview',filesep,'Stage2PostReview',dataSetID,'.mat'],dcDataSetDir)

% Load the data
load([mainDC,filesep,'Stage2ResultsPostReview',filesep,'Stage2PostReview',dataSetID,'.mat'])
load([mainDC,filesep,'Stage2ResultsPreReview',filesep,'Stage2PreReview',dataSetID,'.mat'])

%% Impose neuron rulings and organize data

tic
disp('Organizing the data.')

% Load the first channel to get parameters
[data,~,chanInfo] = load_open_ephys_data_faster([dcDataSetDir,'/100_CH1.continuous']);
sampleRate = chanInfo.header.sampleRate;
nT = length(data);

% Get the number of shanks
nShanks = length(unique(neurList(:,1)));

% Load all the neuron spike times
allSpkData = cell([nShanks,1]);
for iShank = 1:nShanks
    
    % Load the firings
    load([dcDataSetDir,filesep,'rezStage3_shank',num2str(iShank),'.mat']);
    
    % Find the total number of neurons found
    nNeurons = length(unique(rez.st3(:,2)));
    
    % Break out the spike bins for each neuron. 
    spkInds = cell([nNeurons,1]);
    for iNeuron = 1:nNeurons
        spkInds{iNeuron} = rez.st3(rez.st3(:,2) == iNeuron,1)';
    end
    
    % Sometimes no spikes are  placed with a template, resulting in a 
    % neuron with no spikes. If this happens, eliminate it!
    iNeuron = 1;
    while iNeuron <= nNeurons
        if isempty(spkInds{iNeuron})
            spkInds(iNeuron) = [];
            nNeurons = nNeurons - 1;
        else
            iNeuron = iNeuron + 1;
        end
    end
    
    % Break out the spike times for each neuron
    allSpkData{iShank} = cell([nNeurons,1]);
    for iNeuron = 1:nNeurons
        allSpkData{iShank}{iNeuron} = spkInds{iNeuron}/sampleRate;
    end
    
end

% Create a spike time raster of the accepted neurons (only neurons marked
% accepted will be accepted) and include merging
toMerge = zeros([length(mergeOrder),1]);
for iMerge = 1:length(mergeOrder)
    if mergeOrder(iMerge) ~= iMerge
        toMerge(iMerge) = 1;
    end
end
base = find((neurRuling == 1) & (toMerge == 0));
others = find((neurRuling == 1) & (toMerge == 1));
spkData = cell([length(base),1]);
spkDataNeurInd = cell([length(base),1]);
for iBase = 1:length(base)
    spkData{iBase} = allSpkData{neurList(base(iBase),1)}{neurList(base(iBase),2)};
    spkDataNeurInd{iBase} = base(iBase);
end
for iOther = 1:length(others)
    spkData{base == mergeOrder(others(iOther))} = sort([spkData{base == mergeOrder(others(iOther))},allSpkData{neurList(others(iOther),1)}{neurList(others(iOther),2)}]);
    spkDataNeurInd{base == mergeOrder(others(iOther))} = [spkDataNeurInd{base == mergeOrder(others(iOther))},others(iOther)];
end
save([dcDataSetDir,filesep,'spkData.mat'],'spkData')
save([dcDataSetDir,filesep,'spkDataNeurInd.mat'],'spkDataNeurInd')

% Get the tracking and lick data streams
dsFactor = 10; % Set the downsampling factor for these signals
xChan = 3; % x-dimension ADC channel: positive values are toward right sipper (Nick's rig)
yChan = 1; % y-dimension ADC channel: positive values are up on IR camera (Nick's rig)
lChan = 5; % left auditory lick recording (Nick's rig)
rChan = 7; % right auditory lick recording (Nick's rig)
[xdata,timestamps,~] = load_open_ephys_data_faster([dcDataSetDir,'/100_ADC',num2str(xChan),'.continuous']);
nBins = floor(length(xdata)/dsFactor);
timestamps(((nBins*dsFactor) + 1):end) = [];
timestamps = mean(reshape(timestamps,[dsFactor,nBins]));
[ydata,~,~] = load_open_ephys_data_faster([dcDataSetDir,'/100_ADC',num2str(yChan),'.continuous']);
[ldata,~,~] = load_open_ephys_data_faster([dcDataSetDir,'/100_ADC',num2str(lChan),'.continuous']);
[rdata,~,~] = load_open_ephys_data_faster([dcDataSetDir,'/100_ADC',num2str(rChan),'.continuous']);

xdata(((nBins*dsFactor) + 1):end) = [];
xdata = mean(reshape(xdata,[dsFactor,nBins]));
ydata(((nBins*dsFactor) + 1):end) = [];
ydata = mean(reshape(ydata,[dsFactor,nBins]));
ldata(((nBins*dsFactor) + 1):end) = [];
ldata = mean(reshape(ldata,[dsFactor,nBins]));
rdata(((nBins*dsFactor) + 1):end) = [];
rdata = mean(reshape(rdata,[dsFactor,nBins]));

save([dcDataSetDir,filesep,'xy.mat'],'xdata','ydata','timestamps')
save([dcDataSetDir,filesep,'licks.mat'],'ldata','rdata','timestamps')


% Get the med associates events
[maEvents,maTimestamps,~] = load_open_ephys_data_faster([dcDataSetDir,'/all_channels.events']);
save([dcDataSetDir,filesep,'maEvents.mat'],'maEvents','maTimestamps')



% Copy the analysis code into the spike sorting directory
% analysisCode = dir(matDir);
% for iFile = 1:length(analysisCode)
%     if (~strcmp(analysisCode(iFile).name,'.')) && (~strcmp(analysisCode(iFile).name,'..'))
%         copyfile([matDir,filesep,analysisCode(iFile).name],dcDataSetDir)
%     end
% end

% Save the parameters variable
save([dcDataSetDir,filesep,'ksMatlabParams.mat'],'dataSetParams')

disp(['Organizing the data took ',num2str(toc,3),' seconds.'])

%% Transport the files to IU Box

tic
disp('Transfering the results back to IU Box.')

cd(dcDataSetDir)

% Make a list of files to transfer. Send everything that isn't already on
% IU Box.
fileNames = dir(dcDataSetDir);
toSend = ones([length(fileNames),1]);
for iFile = 1:length(fileNames)
    if strcmp(fileNames(iFile).name,'.')
        toSend(iFile) = 0;
    end
    if strcmp(fileNames(iFile).name,'..')
        toSend(iFile) = 0;
    end
    if strcmp(fileNames(iFile).name,'all_channels.events')
        toSend(iFile) = 0;
    end
    if strcmp(fileNames(iFile).name,'Continuous_Data.openephys')
        toSend(iFile) = 0;
    end
    if strcmp(fileNames(iFile).name,'messages.events')
        toSend(iFile) = 0;
    end
    if strcmp(fileNames(iFile).name,'settings.xml')
        toSend(iFile) = 0;
    end
    if length(fileNames(iFile).name) >= 3
        if strcmp(fileNames(iFile).name(1:3),'100')
            toSend(iFile) = 0;
        end
    end
end

% Write the bash file to carry out the transfer
TxtCell = cell(1000,1);

% Fill in the bash script
TxtCell{1,1} = '#!/bin/bash';
TxtCell{2,1} = 'lftp <<SCRIPT';
TxtCell{3,1} = 'set ftps:initial-prot ""';
TxtCell{4,1} = 'set ftp:ssl-force true';
TxtCell{5,1} = 'set ftp:ssl-protect-data true';
TxtCell{6,1} = 'open ftps://ftp.box.com:990';
TxtCell{7,1} = ['user nmtimme@indiana.edu ',IUstring];
TxtCell{8,1} = ['cd "',boxDataSetDir,'"'];
temp = regexp(dcDataSetDir,filesep,'split');
TxtCell{9,1} = ['mkdir "SpikeSorting-',temp{end},'"'];
TxtCell{10,1} = ['cd "SpikeSorting-',temp{end},'"'];
iRow = 11;
for iFile = 1:length(fileNames)
    if toSend(iFile) == 1
        TxtCell{iRow,1} = ['put ',fileNames(iFile).name];
        iRow = iRow + 1;
    end
end
TxtCell{iRow,1} = 'exit';
iRow = iRow + 1;
TxtCell{iRow,1} = 'SCRIPT';
iRow = iRow + 1;
nRows = iRow - 1;
TxtCell(iRow:end,:) = [];

% Name the file transfer list
txtfilename = 'fileTransferListDC2IUBox.txt'; 

% Write the information to the text file
fid = fopen(txtfilename, 'w');
for iRow = 1:size(TxtCell,1)
    fprintf(fid, '%s\r\n', TxtCell{iRow,:});
end
fprintf(fid,'\r\n');
fclose(fid);



% Run the bash script to transfer the files
system(['chmod u+x ',txtfilename]);
system(['dos2unix ',txtfilename]);
% system(['./',txtfilename])
% We had to fix the SSL library this way. If the ftp transfer fails, check
% this to see if it has been updated. See manual for more information.
system(['LD_LIBRARY_PATH=/N/soft/cle5/openssl/1.0.2k/lib:$LD_LIBRARY_PATH ',txtfilename])

disp(['Transfering the data back to IU Box took ',num2str(toc,3),' seconds.'])

end

