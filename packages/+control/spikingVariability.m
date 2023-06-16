function epoch_spike_var = spikingVariability(animal, timeAxis, ...
                                        sessionTypePerBin, spikeRateMatrix)

% this function looks at the variability of neuronal spiking in one animal
% by epoch

% input:
% animal: name of the animal
% timeAxis: midpoints of the time bins
% sessionTypePerBin: the type of session (running/sleeping) each bin is
% spikeRateMatrix: the spike rate matrix generated from getSpikeTrain

% output:
% epoch_spike_var: the variance of the spike rate of all epochs of the
% animal

num_cells = size(spikeRateMatrix,1);
running_spikeTimes =  timeAxis(sessionTypePerBin == 1);
[animal_behavior, throwout_times] = table.behavior.lookup(animal, running_spikeTimes);
%[animal_behavior,unique_times] = behaviors.addBehToTable(animal_behavior);

running_spikerate = spikeRateMatrix(:,sessionTypePerBin == 1);

groups = findgroups(animal_behavior.epoch);
num_groups = numel(unique(groups));

epoch_spike_var = zeros(num_groups, num_cells);

for u = 1:num_groups
    epoch_spike = running_spikerate(:,groups == u);
    epoch_spike_var(u,:) = var(epoch_spike,0,2)';
end

end


