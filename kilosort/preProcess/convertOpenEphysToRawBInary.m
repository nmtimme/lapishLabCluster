function ops = convertOpenEphysToRawBInary(ops)

%fname       = fullfile(ops.root, sprintf('%s.dat', ops.fbinary)); % Nick changed this line to the next one
fname        = ops.fbinary;
fidout      = fopen(fname, 'w');
%
clear fs
load(ops.chanMap,'openEphysChans') % Nick added
for j = 1:ops.Nchan
   %fs{j} = dir(fullfile(ops.root, sprintf('*CH%d_*.continuous', j) )); % Nick changed this line to the next one
   fs{j} = dir(fullfile(ops.root,['100_CH',num2str(openEphysChans(j)),'.continuous']));
end
nblocks = cellfun(@(x) numel(x), fs);
if numel(unique(nblocks))>1
   error('different number of blocks for different channels!') 
end
%
nBlocks     = unique(nblocks);

if ~isempty(ops.refData) % Nick added for referencing
    load(ops.refData)
    if nBlocks ~= 1
        error('We do not know how to handle referencing with multiple blocks.')
    end
end

nSamples    = 1024;  % fixed to 1024 for now!

fid = cell(ops.Nchan, 1);

tic
for k = 1:nBlocks
    for j = 1:ops.Nchan
        fid{j}             = fopen(fullfile(ops.root, fs{j}(k).name));
        % discard header information
        fseek(fid{j}, 1024, 0);
    end
    %
    nsamps = 0;
    flag = 1;
    iRef = 1; % Nick added for referencing
    while 1
        samples = zeros(nSamples * 1000, ops.Nchan, 'int16');
        for j = 1:ops.Nchan
            collectSamps    = zeros(nSamples * 1000, 1, 'int16');
            
            rawData         = fread(fid{j}, 1000 * (nSamples + 6), '1030*int16', 10, 'b');

            nbatches        = ceil(numel(rawData)/(nSamples+6));
            for s = 1:nbatches
                rawSamps = rawData((s-1) * (nSamples + 6) +6+ [1:nSamples]);
                collectSamps((s-1)*nSamples + [1:nSamples]) = rawSamps;
            end
            samples(:,j)         = collectSamps;
        end
                                
        if nbatches<1000
            flag = 0;
        end
        if flag==0
            samples = samples(1:s*nSamples, :);
        end
       
        samples         = samples';
        
        if ~isempty(ops.refData) % Nick added for referencing
            jRef = iRef + size(samples,2) - 1;
            samples = samples - repmat(chanRefint16(iRef:jRef),[size(samples,1),1]);
            iRef = jRef + 1;
        end
        
        if ~isempty(ops.quietChan) % Nick added to allow for quieting channels
            for j = 1:ops.Nchan
                if ops.quietChan(j) == 2
                    samples(j,:) = int16(5*randn(size(samples(j,:))));
                end
            end
        end
        
        if ops.dup2noise == 1 % Nick added to allow noisy channels for shanks with less than 4 good channels
            for j = 2:ops.Nchan
                if isequal(samples(1,:),samples(j,:))
                    samples(j,:) = int16(5*randn(size(samples(j,:))));
                end
            end
        end


        
        fwrite(fidout, samples, 'int16');
        
        nsamps = nsamps + size(samples,2);
        
        if flag==0
            break;
        end
    end
    ops.nSamplesBlocks(k) = nsamps;
    
    for j = 1:ops.Nchan
       fclose(fid{j}); 
    end
    
end
    
fclose(fidout);

toc