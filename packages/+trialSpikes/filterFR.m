function [spikeCountMatrix, spikeRateMatrix, avgFR, areaPerNeuron, cell_index] = ...
         filterFR(spikeCountMatrix, spikeRateMatrix, threshold_fr,...
                  timeAxis, areaPerNeuron, cell_index)


% Filters the spike count matrix based on given firing rate threshold

% get the total number of spikes
total_number_of_spikes = sum(spikeCountMatrix,2);
total_length_of_time = timeAxis(end)-timeAxis(1);

avgFR = total_number_of_spikes/total_length_of_time;
%avgFR = sum(spikeRateMatrix,2)/size(spikeRateMatrix,2);
spikeCountMatrix = spikeCountMatrix(avgFR > threshold_fr,:);
spikeRateMatrix = spikeRateMatrix(avgFR > threshold_fr,:);

areaPerNeuron = areaPerNeuron(avgFR > threshold_fr);
cell_index = cell_index(avgFR > threshold_fr,:);
end

