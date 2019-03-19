function [] = stage1Core(info)
%STAGE1CORE runs the first stage of the kilosort spike sorting on Big Red 
%   2. This stage handles moving data from IU Box to the Data Capacitor 
%   and generating a report on basic features of the data to be reviewed 
%   by the user.
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
%       info.username (string): IU Box username


%% Pull variables out of the info structure

boxDataSetDir = info.boxDataSetDir;
dataSetParams = info.dataSetParams;
dcDataSetDir = info.dcDataSetDir;
dataSetID = info.dataSetID;
IUstring = info.IUstring;
mainDC = info.mainDC;
username = info.username;
suffix = info.suffix;
%% Transfer the data from box to the data capacitor

disp('Starting data transfer from IU Box to the Data Capacitor.')
tic

% Go to the data capacitor folder
matDir = cd(dcDataSetDir);

% Make the bash file to transfer the files
TxtCell = cell(1000,1);

% Fill in the bash script
TxtCell{1,1} = '#!/bin/bash';
TxtCell{2,1} = 'lftp <<SCRIPT';
TxtCell{3,1} = 'set ftps:initial-prot ""';
TxtCell{4,1} = 'set ftp:ssl-force true';
TxtCell{5,1} = 'set ftp:ssl-protect-data true';
TxtCell{6,1} = 'open ftps://ftp.box.com:990';
TxtCell{7,1} = ['user ',username,'@iu.edu ',IUstring];
TxtCell{8,1} = ['cd "',boxDataSetDir,'"'];
TxtCell{9,1} = ['get all_channels' suffix '.events'];
TxtCell{10,1} = ['get Continuous_Data' suffix '.openephys'];
TxtCell{11,1} = ['get messages' suffix '.events'];
TxtCell{12,1} = ['get settings' suffix '.xml'];
iRow = 13;
for i = 1:8
    TxtCell{iRow,1} = ['get 100_ADC',num2str(i),suffix,'.continuous'];
    iRow = iRow + 1;
end
for i = 1:3
    TxtCell{iRow,1} = ['get 100_AUX',num2str(i),suffix,'.continuous'];
    iRow = iRow + 1;
end
for i = 1:dataSetParams{1}
    TxtCell{iRow,1} = ['get 100_CH',num2str(i),suffix,'.continuous'];
    iRow = iRow + 1;
end
TxtCell{iRow,1} = 'exit';
iRow = iRow + 1;
TxtCell{iRow,1} = 'SCRIPT';
iRow = iRow + 1;
nRows = iRow - 1;
TxtCell(iRow:end,:) = [];

% Name the file transfer list
txtfilename = 'fileTransferListIUBox2DC.txt'; 

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


disp(['Finished data transfer from IU Box to the Data Capacitor. It took ',num2str(toc,3),' seconds.'])



%% Generate a report about the data for review by the user

disp('Starting basic data report generation.')
tic

% Set the down sampling factor
nDSFact = 10;

% Set the number of seconds in each data chunk to be recorded
nTChunk = 4;

% Go the Matlab directory
cd(matDir)

% Load each channel
chanInfo = cell([dataSetParams{1},1]);
dataChunk = cell([dataSetParams{1},3]);
power = cell([dataSetParams{1},3]);
dataSD = NaN([dataSetParams{1},1]);
TempMedianRef=cell([1,3]);
FiltdataThr= NaN(dataSetParams{1},3);

for iChan = 1:dataSetParams{1}
    [data,timestamps,chanInfo{iChan}] = load_open_ephys_data_faster([dcDataSetDir,'/100_CH',num2str(iChan),suffix,'.continuous']);
        
    if size(data,2) == 1
        data = data';
    end
    if size(timestamps,2) == 1
        timestamps = timestamps';
    end
    
    % Record the standard deviation of the data
    dataSD(iChan) = std(data);
    
    % Get 4 seconds of data early, middle, and late in the recording. Also
    % get the power spectra for these chunks.
    minT = timestamps(1);
    maxT = timestamps(end);
    tStart = NaN([1,3]);
    tStart(1) = minT + 0.05*(maxT - minT);
    tStart(2) = minT + 0.5*(maxT - minT);
    tStart(3) = minT + 0.95*(maxT - minT);
    for i = 1:3
        tempData = data((timestamps >= tStart(i)) & (timestamps <= (tStart(i) + nTChunk)));
        tempTimestamps = timestamps((timestamps >= tStart(i)) & (timestamps <= (tStart(i) + nTChunk)));
        TempMedianRef{i}(iChan,:)= tempData;
        
        % Down sample the data
        n = length(tempData);
        nSm = floor(n/nDSFact);
        if rem(n,nDSFact) ~= 0
            tempData(((nSm*nDSFact) + 1):end) = [];
            tempTimestamps(((nSm*nDSFact) + 1):end) = [];
        end
        tempData = mean(reshape(tempData,[nDSFact,nSm]));
        tempTimestamps = mean(reshape(tempTimestamps,[nDSFact,nSm]));

        dataChunk{iChan,i} = [tempData;tempTimestamps];
        
        % Calculate the power spectrum
        Fs = 1/(tempTimestamps(2) - tempTimestamps(1));
        T = 1/Fs;
        L = length(tempData);
        Y = fft(tempData);
        P2 = abs(Y/L);
        P1 = P2(1:L/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f = Fs*(0:(L/2))/L;
        power{iChan,i} = [P1;f];
        
        
    end

%     % Down sample the data
%     n = length(data);
%     nSm = floor(n/nDSFact);
%     if rem(n,nDSFact) ~= 0
%         data(((nSm*nDSFact) + 1):end) = [];
%         timestamps(((nSm*nDSFact) + 1):end) = [];
%     end
%     data = mean(reshape(data,[nDSFact,nSm]));
%     timestamps = mean(reshape(timestamps,[nDSFact,nSm]));
%     
%     % Calculate the power spectrum
%     y = fft(data);
%     fs = 1/(timestamps(2) - timestamps(1));
%     n = length(data);
%     f = (0:n-1)*(fs/n);
%     power{iChan} = [abs(y).^2/n;f];
    
       
end

% add a filtered row for each data chunk

% design filter between 300 and 3000Hz
sampleRate = chanInfo{1}.header.sampleRate;
[b,a]=ellip(4,0.1,40,[300 3000]*2/(sampleRate));
for i=1:3
    % get median of channels for referencing the chunks
    MedianRef{i}=median(TempMedianRef{i});
    for iChan=1:dataSetParams{1}
        tempFiltered= filtfilt(b,a,TempMedianRef{i}(iChan,:)-MedianRef{i});
        %downsample
        n = length(tempFiltered);
        nSm = floor(n/nDSFact);
        if rem(n,nDSFact) ~= 0
            tempFiltered(((nSm*nDSFact) + 1):end) = [];
        end
        tempFiltered = mean(reshape(tempFiltered,[nDSFact,nSm]));
        dataChunk{iChan,i}=[dataChunk{iChan,i}; tempFiltered];
        FiltdataThr(iChan,i) = 5 * median(abs(tempFiltered))/0.6745;

    end
end
        
% Save a copy of the stage 1 jobs results in the directory for the data set
% and in the directory with all the results for easy shipment back to local
% machines.
save([dcDataSetDir,filesep,'Stage1PreReview',dataSetID],'chanInfo','dataChunk','power','dataSD')
save([mainDC,filesep,'Stage1ResultsPreReview',filesep,'Stage1PreReview',dataSetID],'chanInfo','dataChunk','power','dataSD','FiltdataThr')

disp(['Finished basic data report generation. It took ',num2str(toc,3),' seconds.'])

%% Direct the user to the report file

disp(['The review file for stage 1 is ',mainDC,filesep,'Stage1ResultsPreReview',filesep,'Stage1PreReview',dataSetID])


end

