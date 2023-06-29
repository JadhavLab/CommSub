function raw = filterFR(raw, threshold_fr)
% [spikeCountMatrix, spikeRateMatrix, avgFR, areaPerNeuron, cell_index] = ...
%         filterFR(spikeCountMatrix, spikeRateMatrix, threshold_fr,...
%                  timeAxis, areaPerNeuron, cell_index)
% Filters the spike count matrix based on given firing rate threshold
% Inputs:
%   raw : struct
%   threshold_fr: threshold firing rate, in Hz
% Outputs:
%   spikeCountMatrix: matrix of spike counts
%   spikeRateMatrix: matrix of spike rates
%   avgFR: average firing rate
%   areaPerNeuron: area of each neuron
%   cell_index: cell index

timeAxis = raw.timeBinMidPoints;

% get the total number of spikes
total_number_of_spikes = sum(raw.spikeCountMatrix,2);
total_length_of_time = timeAxis(end)-timeAxis(1);
raw.avgFR = total_number_of_spikes/total_length_of_time;
%avgFR = sum(spikeRateMatrix,2)/size(spikeRateMatrix,2);

raw.spikeCountMatrix = raw.spikeCountMatrix(raw.avgFR > threshold_fr,:);
raw.spikeRateMatrix  = raw.spikeRateMatrix(raw.avgFR > threshold_fr,:);
raw.areaPerNeuron    = raw.areaPerNeuron(raw.avgFR > threshold_fr);
raw.cell_index       = raw.cell_index(raw.avgFR > threshold_fr,:);
raw.avgFR            = raw.avgFR(raw.avgFR > threshold_fr); % fix June 2023

end
