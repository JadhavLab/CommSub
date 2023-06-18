function [spikeCountMatrix, spikeRateMatrix, avgFR, areaPerNeuron, cell_index] = ...
         filterFR(spikeCountMatrix, spikeRateMatrix, threshold_fr,...
                  timeAxis, areaPerNeuron, cell_index)
% [spikeCountMatrix, spikeRateMatrix, avgFR, areaPerNeuron, cell_index] = ...
%         filterFR(spikeCountMatrix, spikeRateMatrix, threshold_fr,...
%                  timeAxis, areaPerNeuron, cell_index)
% Filters the spike count matrix based on given firing rate threshold
% Inputs:
%   spikeCountMatrix: matrix of spike counts, each row is a neuron
%   spikeRateMatrix: matrix of spike rates, each row is a neuron
%   threshold_fr: threshold firing rate, in Hz
%   timeAxis: time axis, in seconds
%   areaPerNeuron: area of each neuron, brain area string
%   cell_index: cell index, cell id number
% Outputs:
%   spikeCountMatrix: matrix of spike counts
%   spikeRateMatrix: matrix of spike rates
%   avgFR: average firing rate
%   areaPerNeuron: area of each neuron
%   cell_index: cell index


% get the total number of spikes
total_number_of_spikes = sum(spikeCountMatrix,2);
total_length_of_time = timeAxis(end)-timeAxis(1);
avgFR = total_number_of_spikes/total_length_of_time;
%avgFR = sum(spikeRateMatrix,2)/size(spikeRateMatrix,2);

spikeCountMatrix = spikeCountMatrix(avgFR > threshold_fr,:);
spikeRateMatrix  = spikeRateMatrix(avgFR > threshold_fr,:);
areaPerNeuron    = areaPerNeuron(avgFR > threshold_fr);
cell_index       = cell_index(avgFR > threshold_fr,:);
avgFR            = avgFR(avgFR > threshold_fr); % fix June 2023

end

