function [cellOfWindows, cutoff] = jointMake(Hstruct, winSize, varargin)
% Generates a set of windows who are JOINTLY defined by ranges of more than one pattern at a time
% .... creates windows around points that begin to meet a joint quantile threshodd
%
%  This function is spiritually related to generate.mate, except that it uses
%  a less confusing struct to simulateously define everything about the joint
%  patterns and the settings that belond to each pattern
%
% 1st struct level, encodes joint dims
% 2nd struct level are parallel pattern dims

ip = inputParser;
ip.addParameter('patternNames', []); % Provided to help order the cellOfWindows output. Otherwise, if empty, it defaults to the order encountered in the struct.
ip.parse(varargin{:})
Opt = ip.Results;

cellOfWindows = {};
if iscell(winSize)
    if numel(winSize) > 1
        error("Nope")
    end
    winSize = winSize{1};
end

% Levels
if ~isempty(Opt.patternNames)
    paraLevel = Opt.patternNames;
else
    paraLevel = string(fieldnames(Hstruct))';
end
nP         = numel(paraLevel);
jointLevel = string(fieldnames(Hstruct.(paraLevel(1))))';
nT         = numel(jointLevel);

% Aliases
LOW = 1;
HIGH = 2;

% Ensure times on same scale
% jkk
% --------------------------
% Find optimal time resolution
time_details = []; % stores [start stop sampperiod]
p = 0;
for para = paraLevel; p = p + 1;
    j = 0;
    for joint = jointLevel; j = j + 1;
        times = Hstruct.(para).(joint).Htimes;
        time_details = [time_details; times(1) times(end) median(diff(times))];
    end
end
optimal_detail = [min(time_details(:,1)), max(time_details(:,2)), min(time_details(:,3))];
new_time = optimal_detail(1):optimal_detail(3):optimal_detail(2);

% Interpolate event matrices to match this
p = 0;
for para = paraLevel; 
    p = p + 1;
    j = 0;
    for joint = jointLevel; 
        j = j + 1;
        item = Hstruct.(para).(joint);
        times = item.Htimes;
        nTimes = numel(times);
        for field = setdiff(string(fieldnames(item)),"Htimes")'
            if numel(item.(field)) == nTimes
                item.(field) = interp1(times, item.(field), new_time);
            end
        end
        item.Htimes = new_time;
        Hstruct.(para).(joint) = item;
    end
end

% Detect windows
% --------------
p = 0;
cutoffs = [];
for para = paraLevel; 
    p = p + 1;

    series = [];
    j = 0;
    for joint = jointLevel
        j = j + 1;

        item = Hstruct.(para).(joint);

        % Obtain quantile cutoffs
        quantileBin = item.quantileBin;
        if isfield(item, "arrayForQuantile")
            samples_to_quantile = item.arrayForQuantile;
        else
            samples_to_quantile = item.H;
        end

        quantile_bin = item.quantileBin;
        assert(quantileBin(1) < quantileBin(2), "Error: quantileBin(1) must be < quantileBin(2)");
        cutoff(p,  j, :) = quantile(samples_to_quantile, quantileBin);
        cutoffs(p, j, :) = cutoff(p, j,:);

        % Utilize those cutoffs to label matching times
        series(j,:) = item.H > cutoff(p, j, 1) & item.H <= cutoff(p, j, 2);

    end

    % Find the transition points into matching times and trigger windows
    series = sum(series, 1) >= size(series,1); % literally an AND over dimension 1
    if sum(series) == 0
        warning("Pattern " + para + " is empty");
    end
    switches = diff(series);
    switches = switches == 1; % Beginnings of the joint pattern
    times = item.Htimes;
    if iscolumn(times)
        times = times';
    end

    % Setup the official windows around these treiggers
    R = times(switches);
    R = R(:);
    if isscalar(winSize)
        R(:,2) = R(:,1)+winSize;
    elseif numel(winSize) == 2
        R = [R+winSize(1), R+winSize(2)];
    else
        error("winSize either a single number (total window time from the trigger) or two numbers (time in front of and behind the trigger)")
    end

    cellOfWindows{p} = R;

end
