function [spikeRaster] = getNeuronSpike(cell_index, animal,  spikeRateMatrix,  varargin)
% this method takes the name of an animal and the cells that were active
% and make a Nneuron * Timepoints matrix

ip = inputParser;
ip.addParameter('spike_origin', "spikeRateMatrix", @isstring);
ip.parse(varargin{:});
opt = ip.Results;

load(animal + "spikes01.mat");
spikes_data = spikes{1};
nNeuron = size(cell_index,1);
% loop over the cells and sitch the times
spikeRaster = [];

if opt.spike_origin == "spikeRateMatrix"  %the timepoints variable should be present in the workspace
    spikeRaster = spikeRateMatrix; 
else
    for iNeuron = 1:nNeuron
        spike_data = [];
        curr_indices = cell_index(iNeuron,:);
        for epoch = 1:size(spikes_data,2)
            curr_spike = spikes_data{epoch}{curr_indices(1)}{curr_indices(2)};
            spike_data = [spike_data, curr_spike.data(:,6)'];
        end
        spikeRaster = [spikeRaster; spike_data];
    end
end

