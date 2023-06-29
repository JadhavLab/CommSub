function startStop = getEpochTimes(animal)



load(animal + "spikes01.mat");
session = 1;
spikes = spikes{session}; % single day epoches

%%
% find the start and stop spiking time of each epoch (2 * numEpoch)
startStop = zeros(2, numel(spikes));
for epoch = 1:numel(spikes)
    epochSpikes = spikes{epoch};
    start_time = intmax;
    end_time = -1;
    for tetrode = 1:numel(epochSpikes)
        tetrodeEpSpikes = epochSpikes{tetrode};
        for neuron = 1:numel(tetrodeEpSpikes)
            neuronTetEpoch = tetrodeEpSpikes{neuron};
            if ~isempty(neuronTetEpoch) && isfield(neuronTetEpoch,'data')...
                    && ~isempty(neuronTetEpoch.data)
                curr_start = neuronTetEpoch.timerange(1);
                curr_end = neuronTetEpoch.timerange(2);
                
                if curr_start < start_time
                    start_time = curr_start;
                end
                
                if curr_end > end_time
                    end_time = curr_end;
                end
            end
        end
    end
    startStop (1,epoch) = start_time;
    startStop (2,epoch) = end_time;
end
end

