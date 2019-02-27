function [] = stage2Core(info)
%STAGE2CORE runs the second stage of the kilosort spike sorting on Big Red 
%   2. This stage handles preparing data for spike sorting, spike sorting, 
%   and generating a report of the spike sorting results.
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

dataSetParams = info.dataSetParams;
dcDataSetDir = info.dcDataSetDir;
dataSetID = info.dataSetID;
mainDC = info.mainDC;


%% Load the results of the user review of the stage 1 results

% Copy the post user review stage 1 results file into the directory for
% this specific analysis to ease transfer back to IU Box
copyfile([mainDC,filesep,'Stage1ResultsPostReview',filesep,'Stage1PostReview',dataSetID,'.mat'],dcDataSetDir)

% Load the data
load([mainDC,filesep,'Stage1ResultsPostReview',filesep,'Stage1PostReview',dataSetID,'.mat'])


%% Create the mda file

% Load the first channel to get parameters
[data,~,chanInfo] = load_open_ephys_data_faster([dcDataSetDir,'/100_CH1.continuous']);
sampleRate = chanInfo.header.sampleRate;
nT = length(data);
clear data

% Calculate the reference, if desired
if dataSetParams{4} == 1
    % Average of non-noise channels as reference
    tic
    disp('Calculating the reference via average.')
    
    probeChanList = 1:length(probeChanShankID);
    nChannels = length(probeChanList);
    
    % Preallocate space for the data
    chanRef = zeros([1,nT]);
    
    nGoodChannels = 0;
    for iProbeChan = 1:length(probeChanList)
        
        % Find the open ephys channel that corresponds to this probe
        % channel
        oeChan = find(openE2probe == probeChanList(iProbeChan));
        
        if chanRuling(oeChan) == 0
            % Accept the channel
            
            % Load the channel data
            [tempChan,~,~] = load_open_ephys_data_faster([dcDataSetDir,'/100_CH',num2str(oeChan),'.continuous']);
            chanRef = chanRef + tempChan';
            nGoodChannels = nGoodChannels + 1;
        end
    end
    chanRef = chanRef/nGoodChannels;
    chanRefint16 = int16(round(chanRef/chanInfo.header.bitVolts));
    
    % Save the reference data
    save([dcDataSetDir,filesep,'ref_allshanks.mat'],'chanRefint16')
    
    disp(['Calculating the reference took ',num2str(toc,3),' seconds.'])
    
elseif dataSetParams{4} == 2
    % Median of non-noise channels as reference
    tic
    disp('Calculating the reference via median.')
    
    probeChanList = 1:length(probeChanShankID);
    nChannels = length(probeChanList);
    
    % Use the number of channels as the batching factor to ensure adequate
    % memory
    nBatches = nChannels;
    tempSize = ceil(nT/nBatches);
    batchSizes = tempSize*ones([nBatches,1]);
    batchSizes(end) = nT - ((nBatches-1)*tempSize);
    
    % Preallocate space for the data
    chanRef = NaN([1,nT]);
    iT = 1;
    for iBatch = 1:nBatches
        
        % Calculate the stop time index
        jT = iT + batchSizes(iBatch) - 1;
        
        % Preallocate space for this batch
        allData = NaN([nChannels,batchSizes(iBatch)]);
        
        for iProbeChan = 1:length(probeChanList)
            
            % Find the open ephys channel that corresponds to this probe
            % channel
            oeChan = find(openE2probe == probeChanList(iProbeChan));
            
            if chanRuling(oeChan) == 0
                % Accept the channel
                
                % Load the channel data
                [tempData,~,~] = load_open_ephys_data_faster([dcDataSetDir,'/100_CH',num2str(oeChan),'.continuous']);
                allData(iProbeChan,:) = tempData(iT:jT);
            end
        end
        
        % Calculate the median
        chanRef(iT:jT) = median(allData,'omitnan');
        
        % Advance the start time index
        iT = jT + 1;
        
    end
    
    chanRefint16 = int16(round(chanRef/chanInfo.header.bitVolts));
    clear allData
    
    % Save the reference data
    save([dcDataSetDir,filesep,'ref_allshanks.mat'],'chanRefint16')
    
    disp(['Calculating the reference took ',num2str(toc,3),' seconds.'])
    
else
    disp('No reference calculated.')
end




% Create files for each individual shank
waveformWindowSize = NaN([max(probeChanShankID),2]);
vBounds = cell([max(probeChanShankID),1]);
shankGeo = cell([max(probeChanShankID),1]);
shankWaveforms = cell([max(probeChanShankID),1]);
neurFR = cell([max(probeChanShankID),1]);
neurRefVio = cell([max(probeChanShankID),1]);
neurNSpk = cell([max(probeChanShankID),1]);
neurAC = cell([max(probeChanShankID),1]);
neurSTA = cell([max(probeChanShankID),1]);
neurList = [];
neurDynFR = cell([max(probeChanShankID),1]);
raster = cell([max(probeChanShankID),1]);
neurStage1Ruling = cell([max(probeChanShankID),1]);
shankAveWaveforms = cell([max(probeChanShankID),1]);
shankDistWaveforms = cell([max(probeChanShankID),1]);
shankWaveformDifs = cell([max(probeChanShankID),1]);
mergeGuide = cell([max(probeChanShankID),1]);
for iShank = 1:(max(probeChanShankID))
    
    disp(['Starting Shank ',num2str(iShank)])
    
    tic
    % Make a list of the probe channels on this shank
    probeChanList = find(probeChanShankID == iShank);
    tempStage1Ruling = NaN([length(probeChanList),1]);
    for iProbeChan = 1:length(probeChanList)
        
        % Find the open ephys channel that corresponds to this probe
        % channel
        oeChan = find(openE2probe == probeChanList(iProbeChan));
        
        % Record the channel ruling from stage 1 review
        tempStage1Ruling(iProbeChan) = chanRuling(oeChan);
        
        % Report the channel rulings to the user
        if ~isfinite(tempStage1Ruling(iProbeChan))
            error(['Open ephys channel ',num2str(oeChan),' was not reviewed.'])
        elseif tempStage1Ruling(iProbeChan) == 0
            disp(['Open ephys channel ',num2str(oeChan),' was accepted.'])
        elseif tempStage1Ruling(iProbeChan) == 1
            disp(['Open ephys channel ',num2str(oeChan),' was silenced and removed from the spike sorting.'])
        elseif tempStage1Ruling(iProbeChan) == 2
            disp(['Open ephys channel ',num2str(oeChan),' was quieted to low level noise.'])
        end
        
    end
    % Remove silenced channels
    probeChanList(tempStage1Ruling == 1) = [];
    
    % Calculate the number of channels
    nChannels = length(probeChanList);
    
    % If there are no good channels on this shank, skip it
    if nChannels == 0
        disp(['Shank ',num2str(iShank),' had all bad channels, so it is being skipped.'])
    else
        
        % If we had less than 4 good electrodes on this shank, notify the
        % user and add duplicates of the first good electrode to later be
        % replaced with noise.
        smShankFlag = 0;
        smShankGoodChan = [];
        if nChannels <= 3
            disp(['Shank ',num2str(iShank),' had only ',num2str(nChannels),' good channels and at least 4 are required, so low level noise channels will be added.'])
            probeChanList = [probeChanList;probeChanList(1)*ones([4-nChannels,1])];
            smShankGoodChan = [ones([1,nChannels]),zeros([1,4-nChannels])];
            nChannels = length(probeChanList);
            smShankFlag = 1;
        end
        
        % Preallocate space for the data
        geo = NaN([nChannels,2]);
        tempStage1Ruling = NaN([nChannels,1]);
        
        for iProbeChan = 1:length(probeChanList)
            
            % Find the open ephys channel that corresponds to this probe
            % channel
            oeChan = find(openE2probe == probeChanList(iProbeChan));
            
            % Record the geometry for this probe channel
            geo(iProbeChan,:) = probeGeo(probeChanList(iProbeChan),:);
            
            % Record the channel ruling from stage 1 review
            tempStage1Ruling(iProbeChan) = chanRuling(oeChan);
            
        end
        
        % If we had less than 4 good electrodes on this shank, move the
        % noise electrodes far from the good ones
        if smShankFlag == 1
            xNoise = (max(geo(smShankGoodChan == 1,1)) - min(geo(smShankGoodChan == 1,1)))/2 + min(geo(smShankGoodChan == 1,1));
            yNoise = max(geo(smShankGoodChan == 1,2)) + 300;
            for iChan = 1:nChannels
                if smShankGoodChan(iChan) == 0
                    geo(iChan,1) = xNoise;
                    geo(iChan,2) = yNoise;
                    yNoise = yNoise + 100;
                end
            end
        end
        
        % If some channels were quieted to low level noise, make a record
        quietFlag = 0;
        if nnz(tempStage1Ruling == 2) > 0
            quietFlag = 1;
        end
            
        
        
        % Add the stage 1 channel rulings
        neurStage1Ruling{iShank} = tempStage1Ruling;
        
        disp(['Finished organizing shank parameters for shank ',num2str(iShank),'. It took ',num2str(toc,3),' seconds.'])
        
        %%%%%%%
        % Create the config file (this is taken from the Standard_Config_MOVEME.m file)
        %%%%%%%
        
        ops = struct;
        
        % Note: for all binary options, 1=yes, 0=no
        ops.GPU                 = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first) (Chris has 0)
        ops.parfor              = 1; % whether to use parfor to accelerate some parts of the algorithm
        ops.verbose             = 0; % whether to print command line progress
        ops.showfigures         = 0; % whether to plot figures during optimization (Chris has 1)
        
        ops.datatype            = 'openEphys';  % binary ('dat', 'bin') or 'openEphys'
        ops.fbinary             = [dcDataSetDir,filesep,'raw_shank',num2str(iShank),'.dat']; % will be created for 'openEphys'
        ops.fproc               = [dcDataSetDir,filesep,'temp_wh_shank',num2str(iShank),'.dat']; % residual from RAM of preprocessed data
        ops.root                = dcDataSetDir; % 'openEphys' only: where raw files are
        
        if ismember(dataSetParams{4},[1,2]) % Nick added for referencing (see convertOpenEphysToRawBInary.m)
            ops.refData         = [dcDataSetDir,filesep,'ref_allshanks.mat']; % The bits reference time series
        else
            ops.refData         = []; % No referencing
        end
        
        if smShankFlag == 1 % Nick added to replace duplicate electrodes with noise on shanks with less than 4 good electrodes (see convertOpenEphysToRawBInary.m)
            ops.dup2noise = 1; % Replace duplicate electrodes with noise
        else
            ops.dup2noise = [];
        end
        
        if quietFlag == 1 % Nick added to allow for quieting channels instead of silencing them
            ops.quietChan = tempStage1Ruling;
        else
            ops.quietChan = [];
        end
        
        ops.fs                  = sampleRate;        % sampling rate		(omit if already in chanMap file)
        ops.NchanTOT            = nChannels;           % total number of channels (omit if already in chanMap file)
        ops.Nchan               = nChannels;           % number of active channels (omit if already in chanMap file)
        ops.Nfilt               = 32*ceil((dataSetParams{2}*nChannels)/32);           % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32) (Chris set at 32)
        %     ops.Nfilt               = 32;           % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32) (Chris set at 32)
        %     ops.nNeighPC            = min([nChannels,12]); % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12) (Nick: Note, this must be less than or equal to nChannels) (Chris set at 4)
        ops.nNeighPC            = 4; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12) (Nick: Note, this must be less than or equal to nChannels) (Chris set at 4)
        %     ops.nNeigh              = 16; % visualization only (Phy): number of neighboring templates to retain projections of (16) (Chris has 5)
        ops.nNeigh              = 5; % visualization only (Phy): number of neighboring templates to retain projections of (16) (Chris has 5)
        
        % options for channel whitening
        ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)
        ops.nSkipCov            = 1; % compute whitening matrix from every N-th batch (1)
        ops.whiteningRange      = nChannels; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)
        
        % define the channel map as a filename (string) or simply an array
        ops.chanMap             = [dcDataSetDir,filesep,'chanMap_shank',num2str(iShank),'.mat']; % make this file using createChannelMapFile.m
        ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info). Nick: I don't think this is relevant if we sort whole shanks.
        % ops.chanMap = 1:ops.Nchan; % treated as linear probe if a chanMap file
        
        % other options for controlling the model and optimization
        ops.Nrank               = 3;    % matrix rank of spike template model (3)
        ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)
        ops.maxFR               = dataSetParams{5};  % maximum number of spikes to extract per batch (20000)
        ops.fshigh              = 300;   % frequency for high pass filtering
        %     ops.fslow               = 3000;   % frequency for low pass filtering (optional) (Chris left this out)
        ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
        ops.scaleproc           = 200;   % int16 scaling of whitened data
        ops.NT                  = 32*1024+ ops.ntbuff;% this is the batch size (try decreasing if out of memory)
        % for GPU should be multiple of 32 + ntbuff
        
        % the following options can improve/deteriorate results.
        % when multiple values are provided for an option, the first two are beginning and ending anneal values,
        % the third is the value used in the final pass.
        ops.Th               = dataSetParams{6};    % threshold for detecting spikes on template-filtered data ([6 12 12])
        ops.lam              = dataSetParams{7};   % large means amplitudes are forced around the mean ([10 30 30])
        ops.nannealpasses    = 4;            % should be less than nfullpasses (4)
        ops.momentum         = 1./[20 400];  % start with high momentum and anneal (1./[20 1000])
        ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)
        ops.mergeT           = .1;           % upper threshold for merging (.1)
        ops.splitT           = .1;           % lower threshold for splitting (.1)
        
        % options for initializing spikes from data
        ops.initialize      = 'no'; %'fromData' or 'no' (Chris has as 'no')
        ops.spkTh           = dataSetParams{3};      % spike threshold in standard deviations (4)
        ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])
        ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])
        ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])
        ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)
        ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)
        
        % load predefined principal components (visualization only (Phy): used for features)
        dd                  = load([pwd,'/kilosort/configFiles/PCspikes2.mat']); % you might want to recompute this from your own data
        ops.wPCA            = dd.Wi(:,1:7);   % PCs
        
        % options for posthoc merges (under construction)
        ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)
        ops.epu     = Inf;
        
        ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.
        
        
        
        %%%%%%%
        % Create the channel map file (this is taken from createChannelMapFile.m)
        %%%%%%%
        
        connected = true(nChannels, 1);
        chanMap = 1:nChannels;
        openEphysChans = NaN([1,nChannels]); % Nick: I added this because I couldn't figure what channel chanMap does and I didn't want to mess with it
        xcoords = NaN([nChannels,1]);
        ycoords = NaN([nChannels,1]);
        for iProbeChan = 1:nChannels
            openEphysChans(iProbeChan) = find(openE2probe == probeChanList(iProbeChan));
            xcoords(iProbeChan) = probeGeo(probeChanList(iProbeChan),1);
            ycoords(iProbeChan) = probeGeo(probeChanList(iProbeChan),2);
        end
        chanMap0ind = chanMap - 1;
        kcoords   = ones(nChannels,1); % grouping of channels (i.e. tetrode groups) Nick: We're sorting whole shanks.
        
        fs = sampleRate; % sampling frequency
        save([dcDataSetDir,filesep,'chanMap_shank',num2str(iShank),'.mat'], ...
            'chanMap','connected', 'xcoords', 'ycoords', 'kcoords', 'chanMap0ind', 'fs','openEphysChans')
        
        disp(['Finished creating the kilosort parameters and channel mapping files for shank ',num2str(iShank),'.'])
        
        %%%%%%%
        % Run the spike sorting (this is modified from master_file_example_MOVEME.m)
        %%%%%%%
        
        disp(['Starting spike sorting on shank ',num2str(iShank),'.'])
        
        tic; % start timer
        if ops.GPU
            gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
        end
        
        % Loading the data is unnecessary here because it loads in
        % preprocessData
        %     if strcmp(ops.datatype , 'openEphys')
        %         ops = convertOpenEphysToRawBInary(ops);  % convert data, only for OpenEphys
        %     end
        
        % Perform the spike sorting
        save(fullfile(ops.root,['startOps_shank',num2str(iShank),'.mat']), 'ops');
        [rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
        save(fullfile(ops.root,['rezStage1_shank',num2str(iShank),'.mat']), 'rez', 'DATA', 'uproj', '-v7.3');
        rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
        save(fullfile(ops.root,['rezStage2_shank',num2str(iShank),'.mat']), 'rez', '-v7.3');
        rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)
        save(fullfile(ops.root,['rezStage3_shank',num2str(iShank),'.mat']), 'rez', '-v7.3');
        
        % AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this)
        rezAutoMerge = merge_posthoc2(rez);
        save(fullfile(ops.root,['rezAutoMerge_shank',num2str(iShank),'.mat']), 'rezAutoMerge', '-v7.3');
        
        % save python results file for Phy
        rezToPhy(rez, ops.root);
        
        % Rename python results files for Phy to include the shank number
        movefile([dcDataSetDir,filesep,'whitening_mat_inv.npy'],[dcDataSetDir,filesep,'whitening_mat_inv_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'whitening_mat.npy'],[dcDataSetDir,filesep,'whitening_mat_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'templates_ind.npy'],[dcDataSetDir,filesep,'templates_ind_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'templates.npy'],[dcDataSetDir,filesep,'templates_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'template_features.npy'],[dcDataSetDir,filesep,'template_features_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'template_feature_ind.npy'],[dcDataSetDir,filesep,'template_feature_ind_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'spike_times.npy'],[dcDataSetDir,filesep,'spike_times_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'spike_templates.npy'],[dcDataSetDir,filesep,'spike_templates_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'spike_clusters.npy'],[dcDataSetDir,filesep,'spike_clusters_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'similar_templates.npy'],[dcDataSetDir,filesep,'similar_templates_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'pc_features.npy'],[dcDataSetDir,filesep,'pc_features_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'pc_feature_ind.npy'],[dcDataSetDir,filesep,'pc_feature_ind_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'channel_positions.npy'],[dcDataSetDir,filesep,'channel_positions_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'channel_map.npy'],[dcDataSetDir,filesep,'channel_map_shank',num2str(iShank),'.npy']);
        movefile([dcDataSetDir,filesep,'amplitudes.npy'],[dcDataSetDir,filesep,'amplitudes_shank',num2str(iShank),'.npy']);
        
        % remove temporary file
        delete(ops.fproc);
        
        disp(['Finished spike sorting for shank ',num2str(iShank),'. It took ',num2str(toc,3),' seconds.'])
        
        %%%%%%%
        % Generate the report on the clustering
        %%%%%%%
        disp(['Starting report creation for shank ',num2str(iShank),'.'])
        
        % Find the total number of neurons found
        nNeurons = length(unique(rez.st3(:,2)));
        
        % Break out the spike times for each neuron.
        spkInds = cell([nNeurons,1]);
        for iNeuron = 1:nNeurons
            spkInds{iNeuron} = rez.st3(rez.st3(:,2) == iNeuron,1)';
        end
        
        % Sometimes no spikes are placed with a template, resulting in a
        % neuron with no spikes. If this happens, eliminate it!
        iNeuron = 1;
        oldNeurNumber = 1:nNeurons;
        while iNeuron <= nNeurons
            if isempty(spkInds{iNeuron})
                spkInds(iNeuron) = [];
                oldNeurNumber(iNeuron) = [];
                nNeurons = nNeurons - 1;
            else
                iNeuron = iNeuron + 1;
            end
        end
        
        % Get merge suggestion data
        mergeGuide{iShank} = cell([nNeurons,1]);
        for iNeuron = 1:nNeurons
            mergeTo = unique(rezAutoMerge.st3(rez.st3(:,2) == oldNeurNumber(iNeuron),5));
            if isequal(mergeTo,oldNeurNumber(iNeuron))
                mergeTo = [];
            end
            if isequal(mergeTo,0)
                mergeGuide{iShank}{iNeuron} = 0;
            else
                mergeGuide{iShank}{iNeuron} = find(ismember(oldNeurNumber,mergeTo));
            end
        end
        
        % Get example waveforms for each neuron
        tic
        nExWaveforms = 100; % Number of example waveforms
        tStart = 2; % Time window before the spike to plot (in ms)
        tEnd = 4; % Time window after the spike to plot (in ms)
        tDown = 3; % Down sampling factor in spike waveforms (integer great than or equal to 1)
        t1 = 0:(1000/sampleRate):tStart;
        t2 = 0:(1000/sampleRate):tEnd;
        t = [-fliplr(t1(2:end)),t2];
        t = t(1:(tDown*(floor(length(t)/tDown)))); % Trim in preparation for downsampling
        t = mean(reshape(t,[tDown,length(t)/tDown])); % Downsample
        t1 = 0:floor(tStart*sampleRate/1000);
        t2 = 0:floor(tEnd*sampleRate/1000);
        tInd = [-fliplr(t1(2:end)),t2];
        tInd = tInd(1:(tDown*floor(length(tInd)/tDown))); % Trim in preparation for downsampling
        neurWaveforms = cell([nNeurons,1]);
        tempvBounds = NaN([nNeurons,2]);
        randSpk = spkInds;
        for iNeuron = 1:nNeurons
            randSpk{iNeuron}(randSpk{iNeuron} <= -tInd(1)) = []; % Eliminate spikes at the very start of the recording
            randSpk{iNeuron}(randSpk{iNeuron} >= (nT - tInd(end))) = []; % Eliminate spikes at the very end of the recording
            if length(randSpk{iNeuron}) > nExWaveforms
                randSpk{iNeuron} = randSpk{iNeuron}(randperm(length(randSpk{iNeuron}),nExWaveforms));
            end
            neurWaveforms{iNeuron} = NaN([2,floor(length(tInd)/tDown),nChannels,length(randSpk{iNeuron})]); % x/y, time, channels, example waveform
            neurWaveforms{iNeuron}(1,:,:,:) = repmat(t,[1,1,nChannels,length(randSpk{iNeuron})]);
        end
        % Proceed channel by channel to get the waveforms from the raw
        % binary file data
        fidout = fopen(ops.fbinary,'r');
        rawData = fread(fidout,[nChannels,nT],'*int16');
        fclose(fidout);
        for iChannel = 1:nChannels
            
%             % Clarify channel name to be a probe channel number
%             iProbeChan = iChannel;
%             
%             % Find the open ephys channel that corresponds to this probe
%             % channel
%             oeChan = find(openE2probe == probeChanList(iProbeChan));
%             
%             if chanRuling(oeChan) == 0
%                 % Accept the channel
%                 
%                 % Load the channel data into the shank data matrix
%                 [shankData,~,~] = load_open_ephys_data_faster([dcDataSetDir,'/100_CH',num2str(oeChan),'.continuous']);
%                 
%                 % Reorient shankData
%                 shankData = shankData';
%                 
%                 % Subtract reference signal, if necessary
%                 if ismember(dataSetParams{4},[1,2])
%                     shankData = shankData - chanRef;
%                 end
%                 
%             elseif chanRuling(oeChan) == 2
%                 % Quiet the channel
%                 
%                 % Fill the shank data for this channel with small non-zero
%                 % numbers
%                 shankData = (10^(-10))*randn(1,nT);
%                 
%             end
            
            % Convert this channel's raw data to uV
            chanData = chanInfo.header.bitVolts*double(rawData(iChannel,:));
            
            % Add the waveforms for this channel for all neurons and
            % example waveforms
            for iNeuron = 1:nNeurons
                for iExWaveform = 1:length(randSpk{iNeuron})
                    tempData = chanData(tInd + randSpk{iNeuron}(iExWaveform));
                    neurWaveforms{iNeuron}(2,:,iChannel,iExWaveform) = mean(reshape(tempData,[tDown,length(tempData)/tDown]));
                end
            end
        end
        
        % Previous method that required loading all data for the whole
        % shank
%         for iNeuron = 1:nNeurons
%             tempInds = spkInds{iNeuron};
%             tempInds(tempInds <= -tInd(1)) = []; % Eliminate spikes at the very start of the recording
%             if length(tempInds) <= nExWaveforms
%                 randSpk = tempInds;
%             else
%                 randSpk = tempInds(randperm(length(tempInds),nExWaveforms));
%             end
%             tempWaveforms = NaN([2,floor(length(tInd)/tDown),nChannels,length(randSpk)]);
%             for iExWaveform = 1:length(randSpk)
%                 
%                 tempData = shankData(:,tInd + randSpk(iExWaveform));
%                 
%                 for iChannel = 1:nChannels
%                     tempWaveforms(2,:,iChannel,iExWaveform) = mean(reshape(tempData(iChannel,:),[tDown,length(tempData(iChannel,:))/tDown]));
%                 end
%                 
%             end
%             tempWaveforms(1,:,:,:) = repmat(t,[1,1,nChannels,length(randSpk)]);
%             
%             % Put it in the main variable for saving
%             neurWaveforms{iNeuron} = tempWaveforms; % x/y, time, channels, example waveform
%         end
        for iNeuron = 1:nNeurons
            % Find the maximum and minimum voltages
            tempY = neurWaveforms{iNeuron}(2,:,:,:);
            tempvBounds(iNeuron,1) = min(tempY(:));
            tempvBounds(iNeuron,2) = max(tempY(:));
        end
        
        % Calculate average waveforms and difference from individual example
        % waveforms
        tempAveWaveforms = cell([nNeurons,1]);
        tempDistWaveforms = NaN([nNeurons,nChannels,4]);
        for iNeuron = 1:nNeurons
            tempAveWaveforms{iNeuron} = NaN([floor(length(tInd)/tDown),nChannels]);
            for iChannel = 1:nChannels
                tempMean = mean(neurWaveforms{iNeuron}(2,:,iChannel,:),4);
                tempAveWaveforms{iNeuron}(:,iChannel) = tempMean;
                tempDifs = NaN([size(neurWaveforms{iNeuron},4),1]);
                for iExWaveform = 1:size(neurWaveforms{iNeuron},4)
                    tempDifs(iExWaveform) = max(abs(neurWaveforms{iNeuron}(2,:,iChannel,iExWaveform) - tempMean));
                end
                tempDifs = sort(tempDifs);
                tempDistWaveforms(iNeuron,iChannel,1) = mean(tempDifs);
                tempDistWaveforms(iNeuron,iChannel,2) = std(tempDifs);
                tempDistWaveforms(iNeuron,iChannel,3) = tempDifs(ceil(0.5*length(tempDifs)));
                tempDistWaveforms(iNeuron,iChannel,4) = tempDifs(ceil(0.95*length(tempDifs)));
            end
        end
        vBounds{iShank} = tempvBounds;
        shankWaveforms{iShank} = neurWaveforms;
        shankAveWaveforms{iShank} = tempAveWaveforms;
        shankDistWaveforms{iShank} = tempDistWaveforms;
        save([dcDataSetDir,filesep,'waveforms_shank',num2str(iShank),'.mat'],'tempAveWaveforms')
        
        % Cut down the number of example waveforms to make the report smaller
        nWaveformsKeep = 40;
        for iNeuron = 1:nNeurons
            shankWaveforms{iShank}{iNeuron}(:,:,:,(nWaveformsKeep + 1):end) = [];
        end
        
        % Calculate differences between spike waveforms to aid the user in
        % merging
        tempWaveformDifs = NaN([nNeurons,nNeurons,nChannels]);
        for iNeuron = 1:nNeurons
            for jNeuron = setdiff(1:nNeurons,iNeuron)
                for iChannel = 1:nChannels
                    tempWaveformDifs(iNeuron,jNeuron,iChannel) = max(abs(tempAveWaveforms{iNeuron}(:,iChannel) - tempAveWaveforms{jNeuron}(:,iChannel)));
                end
            end
        end
        shankWaveformDifs{iShank} = tempWaveformDifs;
        
        disp(['Finished getting average waveforms for shank ',num2str(iShank),'. It took ',num2str(toc,3),' seconds.'])
        
        % Figure out the spacing and sizing for the waveform plots on the shank
        tic
        overlapFlag = 0;
        iStep = 1;
        stepSize = 0.1; % size of x range increase in units of the geo data.
        yRatio = 9/16; % Ratio of y range to x range
        while overlapFlag == 0
            for iProbeChan = 1:(nChannels - 1)
                for jProbeChan = (iProbeChan + 1):nChannels
                    ix1 = geo(iProbeChan,1) - iStep*stepSize;
                    ix2 = geo(iProbeChan,1) + iStep*stepSize;
                    iy1 = geo(iProbeChan,2) - iStep*stepSize*yRatio;
                    iy2 = geo(iProbeChan,2) + iStep*stepSize*yRatio;
                    jx1 = geo(jProbeChan,1) - iStep*stepSize;
                    jx2 = geo(jProbeChan,1) + iStep*stepSize;
                    jy1 = geo(jProbeChan,2) - iStep*stepSize*yRatio;
                    jy2 = geo(jProbeChan,2) + iStep*stepSize*yRatio;
                    if (jx1 >= ix1) && (jx1 <= ix2) && (jy1 >= iy1) && (jy1 <= iy2) % j's lower left corner is in i's figure
                        overlapFlag = 1;
                    elseif (jx2 >= ix1) && (jx2 <= ix2) && (jy1 >= iy1) && (jy1 <= iy2) % j's lower right corner is in i's figure
                        overlapFlag = 1;
                    elseif (jx1 >= ix1) && (jx1 <= ix2) && (jy2 >= iy1) && (jy2 <= iy2) % j's upper left corner is in i's figure
                        overlapFlag = 1;
                    elseif (jx2 >= ix1) && (jx2 <= ix2) && (jy2 >= iy1) && (jy2 <= iy2) % j's upper right corner is in i's figure
                        overlapFlag = 1;
                    end
                end
            end
            
            if overlapFlag == 0
                iStep = iStep + 1;
            end
        end
        
        % Reduce the maximum size without touching by 5%
        xSize = iStep*stepSize*0.95;
        ySize = iStep*stepSize*yRatio*0.95;
        waveformWindowSize(iShank,:) = [xSize,ySize];
        shankGeo{iShank} = geo;
        
        % Calculate the firing rates of the neurons
        tempFR = NaN([nNeurons,1]);
        for iNeuron = 1:nNeurons
            tempFR(iNeuron) = length(spkInds{iNeuron})/(nT/sampleRate);
        end
        neurFR{iShank} = tempFR;
        
        % Calculate the number of spikes for the neurons
        tempNSpk = NaN([nNeurons,1]);
        for iNeuron = 1:nNeurons
            tempNSpk(iNeuron) = length(spkInds{iNeuron});
        end
        neurNSpk{iShank} = tempNSpk;
        
        % Calculate the number of refractory period violators
        dTThresh = 1; % refractory period violation window in ms
        dTThresh = floor(dTThresh*sampleRate/1000); % Convert to indexes from time
        tempRefVio = NaN([nNeurons,1]);
        for iNeuron = 1:nNeurons
            tempRefVio(iNeuron) = nnz(diff(spkInds{iNeuron}) <= dTThresh);
        end
        neurRefVio{iShank} = tempRefVio;
        
        disp(['Finished site map plot and other basic calculations for shank ',num2str(iShank),'. It took ',num2str(toc,3),' seconds.'])
        
        % Calculate the auto-correlations
        tic
        tWindow = 100; % Maximum time before and after a spike in ms
        binT = 1; % Bin size in ms
        indWindow = floor(tWindow*sampleRate/1000);
        edges = (0:floor(tWindow/binT))*binT;
        edges = [-fliplr(edges),edges(2:end)];
        tempAC = zeros([nNeurons,length(edges) - 1]);
        parfor iNeuron = 1:nNeurons
            tempInds = spkInds{iNeuron};
            for iSpike = 1:length(tempInds)
                smInds = tempInds((tempInds >= (tempInds(iSpike) - indWindow)) & (tempInds <= (tempInds(iSpike) + indWindow))); % Just get the spikes around this main spike
                smInds(smInds == tempInds(iSpike)) = []; % Remove the main spike
                if ~isempty(smInds)
                    smInds = smInds - tempInds(iSpike); % Make indexes relative to the main spike
                    smInds = smInds/(sampleRate/1000); % Convert indexes to time
                    tempAC(iNeuron,:) = tempAC(iNeuron,:) + histcounts(smInds,edges);
                end
            end
        end
        neurAC{iShank} = [edges(1:(end - 1)) + 0.5*(edges(2) - edges(1));tempAC];
        
        disp(['Finished calculating autocorrelations for shank ',num2str(iShank),'. It took ',num2str(toc,3),' seconds.'])
        
        % Calculate the spike triggered average
        tic
        tWindow = 50; % Maximum time before and after a spike in ms
        binT = 1; % Bin size in ms
        indWindow = floor(tWindow*sampleRate/1000);
        edges = (0:floor(tWindow/binT))*binT;
        edges = [-fliplr(edges),edges(2:end)];
        tempAllSTA = zeros([nNeurons,nNeurons,length(edges) - 1]);
        parfor iNeuron = 1:(nNeurons - 1)
            tempiInds = spkInds{iNeuron};
            tempjSTA = zeros([1,nNeurons,length(edges) - 1]);
            for jNeuron = (iNeuron + 1):nNeurons
                tempjInds = spkInds{jNeuron};
                tempSTA = zeros([1,length(edges) - 1]);
                for iSpike = 1:length(tempiInds)
                    smInds = tempjInds((tempjInds >= (tempiInds(iSpike) - indWindow)) & (tempjInds <= (tempiInds(iSpike) + indWindow))); % Just get the spikes around this main spike
                    if ~isempty(smInds)
                        smInds = smInds - tempiInds(iSpike); % Make indexes relative to the main spike
                        smInds = smInds/(sampleRate/1000);
                        tempSTA = tempSTA + histcounts(smInds,edges);
                    end
                end
                tempjSTA(1,jNeuron,:) = tempSTA;
            end
            tempAllSTA(iNeuron,:,:) = tempjSTA;
        end
        for iNeuron = 1:(nNeurons - 1)
            for jNeuron = (iNeuron + 1):nNeurons
                tempAllSTA(jNeuron,iNeuron,:) = flipud(squeeze(tempAllSTA(iNeuron,jNeuron,:)));
            end
        end
        tempAllSTA(1,1,:) = edges(1:(end - 1)) + 0.5*(edges(2) - edges(1)); % Put the time vector where the STA between neuron 1 and 1 would be
        neurSTA{iShank} = tempAllSTA;
        
        disp(['Finished calculating cross-correlations for shank ',num2str(iShank),'. It took ',num2str(toc,3),' seconds.'])
        
        % Note, Nick tried to limit the number of spikes in the spike triggered
        % average calculation, but this produced some odd results and produced
        % estimations that were subject to noise. This was implimented in the
        % code below. The original code is shown above. If there are many
        % spikes (> ~100000), the code below is necessary to allow it to
        % finish.
        %     % Calculate the auto-correlations
        %     tic
        %     tWindow = 100; % Maximum time before and after a spike in ms
        %     binT = 1; % Bin size in ms
        %     nMaxSpikes = 10000; % Set a limit in the number of spikes to examine
        %     indWindow = floor(tWindow*sampleRate/1000);
        %     edges = (0:floor(tWindow/binT))*binT;
        %     edges = [-fliplr(edges),edges(2:end)];
        %     tempAC = zeros([nNeurons,length(edges) - 1]);
        %     parfor iNeuron = 1:nNeurons
        %         tempIndsCenter = spkInds{iNeuron};
        %         tempInds = spkInds{iNeuron};
        %         if length(tempIndsCenter) > nMaxSpikes
        %             tempIndsCenter = tempIndsCenter(randperm(length(tempIndsCenter),nMaxSpikes));
        %         end
        %         for iSpike = 1:length(tempIndsCenter)
        %             smInds = tempInds((tempInds >= (tempIndsCenter(iSpike) - indWindow)) & (tempInds <= (tempIndsCenter(iSpike) + indWindow))); % Just get the spikes around this main spike
        %             smInds(smInds == tempIndsCenter(iSpike)) = []; % Remove the main spike
        %             if ~isempty(smInds)
        %                 smInds = smInds - tempIndsCenter(iSpike); % Make indexes relative to the main spike
        %                 smInds = smInds/(sampleRate/1000); % Convert indexes to time
        %                 tempAC(iNeuron,:) = tempAC(iNeuron,:) + histcounts(smInds,edges);
        %             end
        %         end
        %         if length(tempIndsCenter) < length(tempInds)
        %             tempAC(iNeuron,:) = tempAC(iNeuron,:)*(length(tempInds)/length(tempIndsCenter));
        %         end
        %     end
        %     neurAC{iShank} = [edges(1:(end - 1)) + 0.5*(edges(2) - edges(1));tempAC];
        %
        %     disp(['Finished calculating autocorrelations for shank ',num2str(iShank),'. It took ',num2str(toc,3),' seconds.'])
        %
        %     % Calculate the spike triggered average
        %     tic
        %     tWindow = 50; % Maximum time before and after a spike in ms
        %     binT = 1; % Bin size in ms
        %     nMaxSpikes = 10000; % Set a limit in the number of spikes to examine
        %     indWindow = floor(tWindow*sampleRate/1000);
        %     edges = (0:floor(tWindow/binT))*binT;
        %     edges = [-fliplr(edges),edges(2:end)];
        %     tempAllSTA = zeros([nNeurons,nNeurons,length(edges) - 1]);
        %     parfor iNeuron = 1:(nNeurons - 1)
        %         tempiInds = spkInds{iNeuron};
        %         if length(tempiInds) > nMaxSpikes
        %             tempiInds = tempiInds(randperm(length(tempiInds),nMaxSpikes));
        %         end
        %         tempjSTA = zeros([1,nNeurons,length(edges) - 1]);
        %         for jNeuron = (iNeuron + 1):nNeurons
        %             tempjInds = spkInds{jNeuron};
        %             tempSTA = zeros([1,length(edges) - 1]);
        %             for iSpike = 1:length(tempiInds)
        %                 smInds = tempjInds((tempjInds >= (tempiInds(iSpike) - indWindow)) & (tempjInds <= (tempiInds(iSpike) + indWindow))); % Just get the spikes around this main spike
        %                 if ~isempty(smInds)
        %                     smInds = smInds - tempiInds(iSpike); % Make indexes relative to the main spike
        %                     smInds = smInds/(sampleRate/1000);
        %                     tempSTA = tempSTA + histcounts(smInds,edges);
        %                 end
        %             end
        %             if length(tempiInds) < spkInds{iNeuron}
        %                 tempjSTA(1,jNeuron,:) = tempSTA*(length(spkInds{iNeuron})/length(tempiInds));
        %             else
        %                 tempjSTA(1,jNeuron,:) = tempSTA;
        %             end
        %         end
        %         tempAllSTA(iNeuron,:,:) = tempjSTA;
        %     end
        %     for iNeuron = 1:(nNeurons - 1)
        %         for jNeuron = (iNeuron + 1):nNeurons
        %             tempAllSTA(jNeuron,iNeuron,:) = flipud(squeeze(tempAllSTA(iNeuron,jNeuron,:)));
        %         end
        %     end
        %     tempAllSTA(1,1,:) = edges(1:(end - 1)) + 0.5*(edges(2) - edges(1)); % Put the time vector where the STA between neuron 1 and 1 would be
        %     neurSTA{iShank} = tempAllSTA;
        %
        %     disp(['Finished calculating cross-correlations for shank ',num2str(iShank),'. It took ',num2str(toc,3),' seconds.'])
        
        % Add the neurons for this shank to the list
        tic
        neurList = [neurList;[iShank*ones([nNeurons,1]),(1:nNeurons)']];
        
        % Calculate firing rates throughout recording
        binT = 1000*60; % Bin size in ms
        edges = 0:(binT*(sampleRate/1000)):nT; % The histogram edges
        tempDynFR = NaN([nNeurons,length(edges) - 1]);
        for iNeuron = 1:nNeurons
            tempDynFR(iNeuron,:) = histcounts(spkInds{iNeuron},edges)/(binT/1000);
        end
        edges = edges/(sampleRate*60); % Convert edges values from indexes to seconds
        neurDynFR{iShank} = [edges(1:(end - 1)) + 0.5*(edges(2) - edges(1));tempDynFR];
        
        
        % Get information for plotting raster
        tStart = 10; % Start time of raster graph in seconds
        tEnd = 2500; % End time of raster graph in seconds
        tempRaster = spkInds;
        for iNeuron = 1:nNeurons
            tempRaster{iNeuron} = tempRaster{iNeuron}/sampleRate; % Convert to seconds
            tempRaster{iNeuron}(tempRaster{iNeuron} > tEnd) = [];
            tempRaster{iNeuron}(tempRaster{iNeuron} < tStart) = [];
        end
        raster{iShank} = tempRaster;
        
        
        disp(['Finished firing rate and rast plot calculations for shank ',num2str(iShank),'. It took ',num2str(toc,3),' seconds.'])
        
    end
end

% Save the report for review
disp('Saving report.')

% Save a copy of the stage 2 jobs results in the directory for the data set
% and in the directory with all the results for easy shipment back to local
% machines.
save([dcDataSetDir,filesep,'Stage2PreReview',dataSetID],'waveformWindowSize','vBounds','shankGeo','shankWaveforms','neurFR','neurRefVio','neurNSpk','neurAC','neurSTA','neurList','neurDynFR','raster','neurStage1Ruling','shankAveWaveforms','shankDistWaveforms','shankWaveformDifs','mergeGuide')
save([mainDC,filesep,'Stage2ResultsPreReview',filesep,'Stage2PreReview',dataSetID],'waveformWindowSize','vBounds','shankGeo','shankWaveforms','neurFR','neurRefVio','neurNSpk','neurAC','neurSTA','neurList','neurDynFR','raster','neurStage1Ruling','shankAveWaveforms','shankDistWaveforms','shankWaveformDifs','mergeGuide')


%% Direct the user to the report file

disp(['The review file for stage 2 is ',mainDC,filesep,'Stage2ResultsPreReview',filesep,'Stage2PreReview',dataSetID])






end

