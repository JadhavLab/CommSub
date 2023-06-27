function [spikeSampleMatrix, spikeSampleTensor, trialTimes] = ...
    generate(spikeCountMatrix, timeAxis, cellOfWindows, samplesPerTrial,...
    numResult)
% function [spikeSampleMatrix, spikeSampleTensor] = ...
%     generate(spikeCountMatrix, timeAxis, cellOfWindows, samplesPerTrial,...
%     numResult)
% Generate spike sample matrix and tensor
% 
% Input
% -----
% spikes : numeric, time x neurons
%   spikeCountMatrix or spikeRateMatrix
%
% timeAxis: 1 x  time
%   Times correpsonding to the the time bins of spikes
%
% cellOfWindows : numeric, trials  x 2
%   Start and stop time of each trial
%
% samplesPerTrial : scalar
%   Number of times to sample equidistant per trial

% Output
% ------
% spikeSampleMatrix : {neuron x (time*trial) per network pattern}
% spikeSampleTensor : {neuron x time x trial per network pattern}
%
% Feb 16, 2021 modified to include trialTimes
disp('Generating spike sample matrix and tensor')

[nNeurons,~]      = size(spikeCountMatrix);
spikeSampleTensor = cell(size(cellOfWindows));
spikeSampleMatrix = cell(size(cellOfWindows));
trialTimes        = cell(size(cellOfWindows));
numResult         = numel(cellOfWindows);

for iPattern = progress(1:numResult, 'Title', 'Patterns')

    try
    [nTrials,~] = size(cellOfWindows{iPattern});
    catch
        keyboard
    end
    currTensor = zeros(nNeurons,samplesPerTrial,nTrials);
    
    currPattern   = cellOfWindows{iPattern};
    if isempty(currPattern)
        windows = 0;
    else
        windows = length(currPattern(:,1));
    end
    
    tq = zeros(1,samplesPerTrial);
    trialTimes{iPattern} = nan(windows, samplesPerTrial);
    for i = progress(1:windows, 'Title', 'Times')
        currstart = currPattern(i,1);
        currstop  = currPattern(i,2);
        
        logicalIndexes = timeAxis>currstart & timeAxis<=currstop;
        
        if sum(logicalIndexes) < 2
            continue;
        end
        
        t = timeAxis(logicalIndexes);
        %t_index = find(logicalIndexes);
        
        interval = (currstop-currstart)/samplesPerTrial;
        
        %Generte query points
        tq(1) = currstart+interval/2;
        for n = 2:samplesPerTrial
            tq(n) = tq(n-1) + interval;
        end
        
        % RY : added test assertions
        assert( min(tq) >= currstart, ...
            'violation: query times below window start')
        assert( max(tq) <= currstop,  ...
            'violation: query times above window stop')
        
        % Per neuron interpolate
        for j = 1:nNeurons
            try
                neuron_spiking = spikeCountMatrix(j, logicalIndexes);
            catch
                keyboard
            end
            try
                interped_spiking = interp1(t,neuron_spiking, tq, 'linear');
            catch
                keyboard
                continue
            end
            currTensor(j,:,i) = interped_spiking;
        end

        trialTimes{iPattern}(i,:) = tq;
    end

    spikeSampleTensor{iPattern} = currTensor;
    temp = reshape(currTensor,[nNeurons,samplesPerTrial*nTrials]);
    temp(isnan(temp)) = 0;

    disp("Checking for rank deficiency")
    if rank(temp) < size(temp, 1)
        disp("Rank deficient ... solving ...")
        temp = clean.solveRankDeficient(temp);
        spikeSampleTensor{iPattern} = reshape(temp, ...
            [nNeurons, samplesPerTrial, nTrials]);
    end
        
    spikeSampleMatrix{iPattern} = temp;
end
