function out = triggered_spectrogram(Patterns_overall, Spk, efizz, varargin)
% SPEC_ANALYSIS Get windows of efizz data around comm sub threshold crossings
%   
% Inputs:
%   Patterns_overall: struct array containing the CCA patterns for each
%       comm sub
%   Spk: struct containing the spike data
%   efizz: struct containing the efizz data
%
% Outputs:
%   out: struct containing the average efizz data for each comm sub

ip = inputParser();
ip.addParameter('ploton', false, @islogical);
ip.addParameter('windowsize', 50, @isnumeric);
ip.addParameter('specNames',  {'S1', 'S2', 'Cavg','wpli_avg'}, @iscellstr);
ip.addParameter('components', [1,2,3,4], @isnumeric);
ip.addParameter('figAppend', "", @(x) isstring(x));
ip.addParameter('folder', "triggered_spectrogram", @(x) isstring(x));
ip.addParameter('quantile_threshold', 0.95, @isnumeric);
ip.addParameter('runtype', 1, @isnumeric); % 1 = run, 0 = rest
ip.parse(varargin{:});
Opt = ip.Results;

Opt.folder = string(Opt.folder) + "_run=" + Opt.runtype + filesep;

if ~exist(figuredefine(Opt.folder), 'dir')
    mkdir(figuredefine(Opt.folder));
end
if ~isempty(Opt.folder) && ~startsWith(Opt.folder, filesep)
    Opt.figAppend = "_" + Opt.figAppend;
end
% Choose your quantile threshold
quantile_threshold = Opt.quantile_threshold;
% Define the window size around the threshold crossing (e.g., 10 time bins
% before and after)
window_size = Opt.windowsize;
% Determine total window size
total_window_size = 2 * window_size + 1;  % accounting for the point itself
% Find indices of pattern with name="overall"
% overall_ind = find(isequal([Patterns_overall.name], "Overall"));

disp("Starting with " + length(Patterns_overall) + " patterns and " + ...
    length(Opt.components) + " components and ploton = " + Opt.ploton);

Opt.means = [];
Opt.mins = [];
for field = Opt.specNames
    if contains(field, "wpli")
        efizz.(field{:})(isnan(efizz.(field{:}))) = 0;
    end
    Opt.means.(field{:}) = mean(efizz.(field{:}), 1);
    Opt.mins.(field{:}) = min(efizz.(field{:}), [], 1);
end

% Loop over all patterns
for i = progress(1:numel(Patterns_overall), 'Title', 'Patterns')

    [d, p] = ind2sub(size(Patterns_overall), i);

    if isempty(Patterns_overall(i).cca) || isempty(Patterns_overall(i).cca.u) || ...
       isempty(Patterns_overall(i).cca.v) || ...
        Patterns_overall(i).directionality == "hpc-hpc" 
        continue
    end

    % Extract u and v from the current pattern
    runs = Spk.sessionTypePerBin==Opt.runtype;
    a = Patterns_overall(i).cca.a;
    b = Patterns_overall(i).cca.b;
    ca1 = Spk.areaPerNeuron == "CA1";
    pfc = Spk.areaPerNeuron == "PFC";
    source = Spk.spikeRateMatrix(ca1, runs);
    target = Spk.spikeRateMatrix(pfc, runs);
    u = source' * a;
    v = target' * b;
    sp_time = Spk.timeBinMidPoints(runs);
    assert(numel(sp_time) == size(u,1), ... 
        'Spike time and u must be the same length');

    % Calculate the quantile threshold for u and v
    u_threshold = quantile(u, quantile_threshold);
    v_threshold = quantile(v, quantile_threshold);

    % Find the time bins where u or v cross their respective thresholds
    u_above_threshold = u > u_threshold;
    v_above_threshold = v > v_threshold;

    if isempty(Patterns_overall(i).nameFull)
        name = Patterns_overall(i).name;
    else
        name = Patterns_overall(i).nameFull;
    end
    all_threshold_crossed = u_above_threshold & v_above_threshold;

    for comp = progress(Opt.components(:)','Title','Components')

        threshold_crossed = all_threshold_crossed(:,comp);  % make sure it's a column vector

        % Determine the corresponding times in efizz
        threshold_crossed_times = sp_time(threshold_crossed); 

        % Find the indices of these times in the efizz data
        efizz_indices = interp1(efizz.t, 1:length(efizz.t), threshold_crossed_times, 'nearest');
        spike_indices = find(threshold_crossed);
        efizz_indices = efizz_indices(:);
        spike_indices = spike_indices(:);
        both = [efizz_indices, spike_indices];
        both = both(~isnan(both(:,1)),:);
        efizz_indices = both(:,1);
        % spike_indices = both(:,2);
        % Find start/stop times for spike windows
        inds_starts = efizz_indices-window_size;
        inds_stops = efizz_indices+window_size;
        % Throw away out of bounds indices
        both = both(inds_starts > 0 & inds_stops <= length(efizz.t),:);
        efizz_indices = both(:,1);
        spike_indices = both(:,2);
        % Find start/stop times for spike windows
        sp_start_indices = interp1(sp_time, 1:numel(sp_time), efizz.t(efizz_indices+window_size), 'nearest');
        sp_start_indices = sp_start_indices(:);
        sp_stop_indices  = interp1(sp_time, 1:numel(sp_time), efizz.t(efizz_indices+window_size), 'nearest');
        sp_stop_indices  = sp_stop_indices(:);
        sp_win_size      = mode(abs([sp_stop_indices - spike_indices; spike_indices - sp_start_indices]))+1;
        total_sp_win_size = 2 * sp_win_size + 1;

        % Initialize empty matrices to store the data segments for each efizz
        % variable
        freqsize = size(efizz.S1, 2);
        segments = cell(1, length(Opt.specNames));
        for j = 1:length(Opt.specNames)
            segments{j} = nan(total_window_size, freqsize, length(efizz_indices));
        end
        time_segments    = nan(total_window_size, length(efizz_indices));
        u_segments       = nan(total_sp_win_size, size(u, 2), length(efizz_indices));
        v_segments       = nan(total_sp_win_size, size(v, 2), length(efizz_indices));
        sp_time_segments = nan(total_sp_win_size, length(efizz_indices));

        % For each efizz index, grab the data in a window around that index
        % for j = progress(1:min(size(both,1),10),'Title','Extracting efizz data')
        for j = progress(1:size(both,1),'Title','Extracting efizz data')
            index = efizz_indices(j);
            spike_index = spike_indices(j);
            eftimes = efizz.t(index-window_size:index+window_size);
            if any(diff(eftimes) > 2)
                continue
            end
            if index-window_size > 0 && index+window_size <= size(efizz.S1, 1)  ...
                && spike_index-sp_win_size > 0 && spike_index+sp_win_size <= size(u, 1)
                % Make sure the window does not exceed the matrix dimensions
                for k = 1:length(Opt.specNames)
                    segments{k}(:,:,j) = ...
                        efizz.(Opt.specNames{k})(index-window_size:index+window_size, :);
                end
                time_segments(:,j) = linspace(efizz.t(index-window_size), ...
                                              efizz.t(index+window_size), ...
                                              total_window_size);
                % Spiking
                sp_time_segments(:,j) = ...
                    sp_time(spike_index-sp_win_size:spike_index+sp_win_size);
                u_segments(:,:,j) = ...
                    u(spike_index-sp_win_size:spike_index+sp_win_size, :);
                v_segments(:,:,j) = ...
                    v(spike_index-sp_win_size:spike_index+sp_win_size, :);
            end
        end

        disp("Averaging " + length(efizz_indices) + " segments for " + name + ...
             " direction " + d + " comp " + comp);
        for j = 1:length(Opt.specNames)
            s=segments{j};
            s=s(:,:,~isnan(s(1,1,:)));
            spec_avg.(Opt.specNames{j}) = mean(s, 3);
        end
        time_segments = time_segments(:,~isnan(time_segments(1,:)));
        time_avg = mean(time_segments, 2);

        % Length
        u_segments = u_segments(:,:,~isnan(u_segments(1,1,:)));
        v_segments = v_segments(:,:,~isnan(v_segments(1,1,:)));
        u_average  = mean(u_segments, 3);
        v_average  = mean(v_segments, 3);

        out(i,comp).spec_avg = spec_avg;
        out(i,comp).threshold_crossed_times = threshold_crossed_times;
        out(i,comp).u_average = u_average;
        out(i,comp).v_average = v_average;
        out(i,comp).u_threshold = u_threshold(comp);
        out(i,comp).v_threshold = v_threshold(comp);
        out(i,comp).time_avg = time_avg;
        out(i,comp).name = name;
        out(i,comp).comp = comp;
        out(i,comp).direction = d;

    end

    if Opt.ploton
        disp("Plotting " + name + " direction " + d);
        plots.triggered_spec_struct(out(i,:), efizz, Opt);
    end
end
