function out = event_analysis(Patterns_overall, Spk, Events)
%EVENT_ANALYSIS   Calculate CCA r-values for each event in a given pattern.
%
%  out = event_analysis(Components, Patterns_overall, Spk, Events)
%
%  INPUTS:
%  Patterns_overall
%  Spk:  Spike data structure 
%  Events:  Structure containing event times
%

% Assuming 'area1' and 'area2' are the indices of the areas you're interested in
area1 = find(strcmp(Spk.areaPerNeuron, 'PFC'));
area2 = find(strcmp(Spk.areaPerNeuron, 'HPC'));

% Create an empty matrix to store the CCA r-values for each event
event_r_values = cell(length(Events.times), 1);
event_u = cell(length(Events.times), 1);
event_v = cell(length(Events.times), 1);
out = struct(...
    'event_r_values', event_r_values, ...
    'event_u', event_u, ...
    'event_v', event_v);
out = repmat(out, size(Patterns_overall));

% Loop over all patterns
for i = 1:length(Patterns_overall)

    % Extract a and b from the current pattern
    a = Patterns_overall(i).cca.a;
    b = Patterns_overall(i).cca.b;

    % Loop over all events
    for j = 1:length(Events.times)

        % Find the time bins that correspond to the current event
        event_time_bins = find(Spk.timeBinStartEnd >= Events.times(j,1) & Spk.timeBinStartEnd <= Events.times(j,2));

        % Extract the spike data for the two areas during this event
        area1_spikes = Spk.spikeCountMatrix(area1, event_time_bins);
        area2_spikes = Spk.spikeCountMatrix(area2, event_time_bins);

        % Project the spike data onto the space defined by a and b to get u and v
        u = area1_spikes' * a;
        v = area2_spikes' * b;

        % Calculate the correlation between u and v
        r = corr(u, v);
        
        % Store the CCA r-value and canonical variates for this event
        event_r_values(j) = r;  % assuming we are interested in the first canonical correlation
        event_u{j} = u;
        event_v{j} = v;
    end

    % Store the CCA r-values and canonical variates for this pattern
    out(i).event_r_values = event_r_values;
    out(i).event_u        = event_u;
    out(i).event_v        = event_v;

end

