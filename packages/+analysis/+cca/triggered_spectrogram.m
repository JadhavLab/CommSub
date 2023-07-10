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
ip.addParameter('specNames', {'S1', 'S2', 'Cavg','Ctoppair'}, @iscellstr);
ip.addParameter('components', [1,2,3,4,5,6,7,8], @isnumeric);
ip.addParameter('figAppend', "", @(x) isstring(x));
ip.addParameter('folder', "triggered_spectrogram", @(x) isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;

if ~exist(figuredefine(Opt.folder), 'dir')
    mkdir(figuredefine(Opt.folder));
end

% if isempty(Patterns_overall) || isempty(Patterns_overall(1).cca) || ...
%    isempty(Patterns_overall(1).cca.u) || isempty(Patterns_overall(1).cca.v)
%     out = [];
%     return
% end

% Choose your quantile threshold
quantile_threshold = 0.95;

% Define the window size around the threshold crossing (e.g., 10 time bins
% before and after)
window_size = Opt.windowsize;

% Determine total window size
total_window_size = 2*window_size + 1;  % accounting for the point itself

% Loop over all patterns
for i = progress(1:length(Patterns_overall), 'Title', 'Patterns')

    [d, n] = ind2sub(size(Patterns_overall), i);

    if isempty(Patterns_overall(i).cca) || isempty(Patterns_overall(i).cca.u) || ...
       isempty(Patterns_overall(i).cca.v)
        continue
    end

    % Extract u and v from the current pattern
    u = Patterns_overall(i).cca.u;
    v = Patterns_overall(i).cca.v;

    % Calculate the quantile threshold for u and v
    u_threshold = quantile(u, quantile_threshold);
    v_threshold = quantile(v, quantile_threshold);

    % Find the time bins where u or v cross their respective thresholds
    u_above_threshold = u > u_threshold;
    v_above_threshold = v > v_threshold;

    % Combine the logical arrays to find when either u or v cross their
    % threshold
    threshold_crossed = u_above_threshold & v_above_threshold;
    for comp = progress(Opt.components(:)','Title','Components')
        threshold_crossed = threshold_crossed(:,comp);  % make sure it's a column vector

        % Determine the corresponding times in efizz
        threshold_crossed_times = Spk.timeBinMidPoints(threshold_crossed); 

        % Find the indices of these times in the efizz data
        efizz_indices = interp1(efizz.t, 1:length(efizz.t), threshold_crossed_times, 'nearest');
        spike_indices = find(threshold_crossed);

        % Initialize empty matrices to store the data segments for each efizz
        % variable
        segments = cell(1, length(Opt.specNames));
        u_segments    = [];
        v_segments    = [];

        % For each efizz index, grab the data in a window around that index
        for j = progress(1:length(efizz_indices),'Title','Extracting efizz data')
            index = efizz_indices(j);
            spike_index = spike_indices(j);
            if index-window_size > 0 && index+window_size <= size(efizz.S1, 1)  % Make sure the window does not exceed the matrix dimensions
                for k = 1:length(Opt.specNames)
                    segments{k} = cat(3, segments{k}, efizz.(Opt.specNames{k})(index-window_size:index+window_size, :));
                end
                u_segments    = cat(3, u_segments, u(spike_index-window_size:spike_index+window_size, :));
                v_segments    = cat(3, v_segments, v(spike_index-window_size:spike_index+window_size, :));
            end
        end

        for j = 1:length(Opt.specNames)
            spec_avg.(Opt.specNames{j}) = mean(segments{j}, 3);
        end
        u_average    = mean(u_segments, 3);
        v_average    = mean(v_segments, 3);

        % Now you have matrices where the third dimension represents different
        % instances when the threshold was crossed. You can now take an average
        % over the third dimension
        if Opt.ploton
            thresholds=struct('u', u_threshold(Opt.comp), 'v', v_threshold(Opt.comp));
            uv = struct('u', u_average(:,Opt.comp), 'v', v_average(:,Opt.comp));
            fig("spec avg direction " + d + " comp " + n); clf;
            plots.triggered_spectrogram(efizz, spec_avg, uv, 'thresholds', thresholds,...
                'nolog', ["Cavg", "Ctoppair"]);
            sz = get(0, 'Screensize');
            pos = mod(Opt.comp, 10) * 0.1;
            set(gcf, 'Position', [pos*sz(3) 0*sz(4) 0.1*sz(3) 1*sz(4)]);
        end

        out(i,comp).S1 = S1_average;
        out(i,comp).S2 = S2_average;
        out(i,comp).C = C_average;
        out(i,comp).wpli = wpli_average;
        out(i,comp).threshold_crossed_times = threshold_crossed_times;
        out(i,comp).u_thresh = u_threshold(Opt.comp);
        out(i,comp).v_thresh = v_threshold(Opt.comp);
    end
end
