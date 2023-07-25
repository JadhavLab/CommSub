function out = event_analysis(Patterns_overall, Spk, Events, Option, varargin)
%EVENT_ANALYSIS   Calculate CCA r-values for each event in a given pattern.
%
%  out = event_analysis(Components, Patterns_overall, Spk, Events)
%
%  INPUTS:
%  Patterns_overall
%  Spk:  Spike data structure 
%  Events:  Structure containing event times
%  Option:  Structure containing options
%   varargin:  'precompute', 0 or 1 (default 0)
%              'method', 'zscore' or 'spikerate' (default 'zscore')
%
%  OUTPUTS:
%  out:  Structure containing CCA r-values for each event in each pattern

ip = inputParser();
ip.addParameter('precompute', 0, @isscalar);
ip.addParameter('method', 'zscore', @ischar);
ip.parse(varargin{:});
Opt = ip.Results;

% Assuming 'area1' and 'area2' are the indices of the areas you're interested in
area1 = find(strcmp(Spk.areaPerNeuron, 'CA1'));
area2 = find(strcmp(Spk.areaPerNeuron, 'PFC'));

% Create an empty matrix to store the CCA r-values for each event
out = struct('W', []);
% out = repmat(out, size(Patterns_overall));

szPatterns = size(Patterns_overall);

% Loop over all patterns
for i = progress(1:numel(Patterns_overall), 'Title', 'Event analysis')

    event_u = cell(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    event_v = cell(length(Events.cellOfWindows), length(Events.cellOfWindows{1}));
    event_r_values = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}), 3);
    event_u_values = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}), 3);
    event_v_values = nan(length(Events.cellOfWindows), length(Events.cellOfWindows{1}), 3);

    if isempty(Patterns_overall(i).cca)
        continue;
    end

    % Extract a and b from the current pattern
    if isfield(Patterns_overall(i).cca, 'a') && Opt.precompute >= 1
        disp("Using precomputed...")
        a = Patterns_overall(i).cca.a;
        b = Patterns_overall(i).cca.b;
    else % if not precomputed
        disp("Computing CCA")
        source = Patterns_overall(i).X_source;
        target = Patterns_overall(i).X_target;
        index_source = Patterns_overall(i).index_source;
        index_target = Patterns_overall(i).index_target;
        if Opt.method == "zscore" && ~munge.detectZscore(source)
            % zscore the fioring
            fprintf("...zscore\n")
            source = zscore(source, 0, 2);
            target = zscore(target, 0, 2);
        elseif Opt.method == "spikerate" && munge.detectZscore(source)
            fprintf("...un-zscore\n")
            % undo z-score normalization             
            muFR   = Spk.muFR(Spk.hpc.to_original(index_source));
            stdFR  = Spk.stdFR(Spk.hpc.to_original(index_source));
            source = source .* stdFR + muFR;
            muFR   = Spk.muFR(Spk.pfc.to_original(index_target));
            stdFR  = Spk.stdFR(Spk.pfc.to_original(index_target));
            target = target .* stdFR + muFR;
        end
        [a,b] = canoncorr(source', target');
    end

    [i1,i2] = ind2sub(szPatterns, i);
    i = {i1, i2};
    if i2 <= Option.nPatternAndControl
        p = i2;
    else
        p = 1:Option.nPatternAndControl;
    end
    directionality = Patterns_overall(i{:}).directionality;

    % Loop over all events
    W = struct();
    for w = p
        windows = Events.cellOfWindows{w};
        for j = 1:length(windows)

            % Find the time bins that correspond to the current event
            event_time_bins = find(Spk.timeBinMidPoints >= windows(j,1) &...
                Spk.timeBinMidPoints <= windows(j,2));
            if isempty(event_time_bins)
                continue;
            end

            % Extract the spike data for the two areas during this event
            if directionality == "hpc-pfc"
                area1_spikes = Spk.spikeCountMatrix(area1, event_time_bins);
                area2_spikes = Spk.spikeCountMatrix(area2, event_time_bins);
            elseif directionality == "hpc-hpc"
                continue;
            end

            % Project the spike data onto the space defined by a and b to get u and v
            u = area1_spikes' * a(:,1:3);
            v = area2_spikes' * b(:,1:3);

            % Calculate the correlation between u and v
            r = [corr(u(:,1), v(:,1)),...
                 corr(u(:,2), v(:,2)),...
                 corr(u(:,3), v(:,3))];
            
            % Store the CCA r-value and canonical variates for this event
            event_r_values(w,j,:) = r;  % assuming we are interested in the first canonical correlation
            event_u_values(w,j,:) = mean(u,1);
            event_v_values(w,j,:) = mean(v,1);
            event_u{w,j} = u;
            event_v{w,j} = v;
        end
    end

    % Store the CCA r-values and canonical variates for this pattern
    out(i{:}).event_r_values = event_r_values;
    out(i{:}).event_u_values = event_u_values;
    out(i{:}).event_v_values = event_v_values;
    out(i{:}).event_u        = event_u;
    out(i{:}).event_v        = event_v;

end
