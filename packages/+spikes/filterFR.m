function [spikeCountMatrix_filtered, spikeRateMatrix_filtered, avgFR] = ...
         filterFR(spikeCountMatrix, spikeRateMatrix, threshold_fr)
% Filters the spike count matrix based on given firing rate threshold
% also returns the firing rate of each neuron so as to keep track 

avgFR = sum(spikeRateMatrix,2)/size(spikeRateMatrix,2);

spikeCountMatrix_filtered = spikeCountMatrix(avgFR > threshold_fr,:);
spikeRateMatrix_filtered  = spikeRateMatrix(avgFR  > threshold_fr,:);
avgFR = avgFR(avgFR > threshold_fr,:);
