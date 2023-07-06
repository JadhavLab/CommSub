function out = ...
    spec_analysis(Patterns_overall, Spk, efizz)
%SPEC_ANALYSIS Get windows of efizz data around comm sub threshold crossings
%   
% Inputs:
%   Patterns_overall: struct array containing the CCA patterns for each
%       comm sub
%   Spk: struct containing the spike data
%   efizz: struct containing the efizz data
%
% Outputs:
%   out: struct containing the average efizz data for each comm sub

% Choose your quantile threshold
quantile_threshold = 0.95;

% Define the window size around the threshold crossing (e.g., 10 time bins before and after)
window_size = 10;

% Determine total window size
total_window_size = 2*window_size + 1;  % accounting for the point itself

% Loop over all patterns
for i = 1:length(Patterns_overall)

    % Extract u and v from the current pattern
    u = Patterns_overall(i).cca.u;
    v = Patterns_overall(i).cca.v;

    % Calculate the quantile threshold for u and v
    u_threshold = quantile(u, quantile_threshold);
    v_threshold = quantile(v, quantile_threshold);

    % Find the time bins where u or v cross their respective thresholds
    u_above_threshold = u > u_threshold;
    v_above_threshold = v > v_threshold;

    % Combine the logical arrays to find when either u or v cross their threshold
    threshold_crossed = u_above_threshold | v_above_threshold;

    % Determine the corresponding times in efizz
    threshold_crossed_times = Spk.timeBinMidPoints(threshold_crossed); 

    % Find the indices of these times in the efizz data
    [~, efizz_indices] = ismember(threshold_crossed_times, efizz.time);

    % Initialize empty matrices to store the data segments for each efizz variable
    S1_segments = [];
    S2_segments = [];
    C_segments = [];
    wpli_segments = [];

    % For each efizz index, grab the data in a window around that index
    for j = 1:length(efizz_indices)
        index = efizz_indices(j);
        if index-window_size > 0 && index+window_size <= size(efizz.S1, 1)  % Make sure the window does not exceed the matrix dimensions
            S1_segments = cat(3, S1_segments, efizz.S1(index-window_size:index+window_size, :));
            S2_segments = cat(3, S2_segments, efizz.S2(index-window_size:index+window_size, :));
            C_segments = cat(3, C_segments, efizz.C(index-window_size:index+window_size, :));
            wpli_segments = cat(3, wpli_segments, efizz.wpli(index-window_size:index+window_size, :));
        end
    end

    % Now you have matrices where the third dimension represents different instances when the threshold was crossed. You can now take an average over the third dimension

    S1_average = mean(S1_segments, 3);
    S2_average = mean(S2_segments, 3);
    C_average = mean(C_segments, 3);
    wpli_average = mean(wpli_segments, 3);

    out(i).S1 = S1_average;
    out(i).S2 = S2_average;
    out(i).C = C_average;
    out(i).wpli = wpli_average;
end

