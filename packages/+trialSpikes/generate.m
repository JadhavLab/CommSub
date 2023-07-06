function Spk = generate(Spk, Events, Option)
% function out = generate(Spk, Option)
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

spikeCountMatrix = Spk.spikeCountMatrix;
timeAxis         = Spk.timeBinMidPoints;
cellOfWindows    = Events.cellOfWindows;
samplesPerTrial  = cast(Option.timesPerTrial, class(timeAxis));
numResult        = Option.nPatternAndControl;

disp('Generating spike sample matrix and tensor')

[nNeurons,~]      = size(spikeCountMatrix);
spikeSampleTensor = cell(size(cellOfWindows));
spikeSampleMatrix = cell(size(cellOfWindows));
trialTimes        = cell(size(cellOfWindows));
numResult         = numel(cellOfWindows);
rowsRemoved       = cell(size(cellOfWindows));
good_neurons       = cell(size(cellOfWindows));

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
        [~, good_neurons{iPattern}] = clean.solveRankDeficient(temp);
    else
        good_neurons{iPattern} = 1:size(temp, 1);
    end
        
    spikeSampleMatrix{iPattern} = temp;

end

% Unify good neurons across patterns
gn = good_neurons{1};
for iPattern = 2:numResult
    gn = intersect(gn, good_neurons{iPattern});
end
good_neurons = gn;

Spk.spikeSampleMatrix = spikeSampleMatrix;
Spk.spikeSampleTensor = spikeSampleTensor;
Spk.trialTimes        = trialTimes;
Spk.rowsRemoved       = rowsRemoved;

% Remove bad neurons from existing variables
Spk.good_neurons_trialPhase      = good_neurons;
Spk.hpc.X = {}; % spikeSampleMatrix
Spk.hpc.T = {}; % spikeSampleTensor
Spk.pfc.X = {};
Spk.pfc.T = {};
for iPattern = 1:numResult
    % Assign spikeSampleMatrix and spikeSampleTensor neurons (rows)
    % to hpc and pfc variables accounting for bad neurons
    % per pattern, select proper areas using Spk.areaPerNeuron 
    % and good_neurons{iPattern}
    % hpc
    hpc_neurons = Spk.areaPerNeuron == "CA1";
    hpc_neurons = intersect(find(hpc_neurons),good_neurons);
    Spk.hpc.X{iPattern} = Spk.spikeSampleMatrix{iPattern}(hpc_neurons, :);
    Spk.hpc.T{iPattern} = Spk.spikeSampleTensor{iPattern}(hpc_neurons, :, :);
    % pfc
    pfc_neurons = Spk.areaPerNeuron == "PFC";
    pfc_neurons = intersect(find(pfc_neurons),good_neurons);
    Spk.pfc.X{iPattern} = Spk.spikeSampleMatrix{iPattern}(pfc_neurons, :);
    Spk.pfc.T{iPattern} = Spk.spikeSampleTensor{iPattern}(pfc_neurons, :, :);
end
Spk.hpc.FR = Spk.avgFR(hpc_neurons);
Spk.pfc.FR = Spk.avgFR(pfc_neurons);
Spk.cell_index = Spk.cell_index(good_neurons,:);
Spk.avgFR = Spk.avgFR(good_neurons);
for iPattern = 1:numResult
    Spk.spikeSampleMatrix{iPattern} = Spk.spikeSampleMatrix{iPattern}(good_neurons,:);
    Spk.spikeSampleTensor{iPattern} = Spk.spikeSampleTensor{iPattern}(good_neurons,:,:);
end
Spk.spikeCountMatrix = Spk.spikeCountMatrix(good_neurons,:);
Spk.spikeRateMatrix  = Spk.spikeRateMatrix(good_neurons,:);
Spk.areaPerNeuron    = Spk.areaPerNeuron(good_neurons);


end
