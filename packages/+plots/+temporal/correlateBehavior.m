function correlateBehavior(Components, Option, behavior, varargin)
% correlateBehavior(Components, Behavior, varargin)
%
%   Correlate the behavior with the components of the model.
%
%   Inputs:

ip = inputParser();
ip.addParameter('behavior', []);
ip.addParameter('unique_times', []);
ip.addParameter('throwout_times', []);
ip.addParameter('use', 'smooth');
ip.addParameter('names', []);
ip.addParameter('figAppend', "", @(x) isstring(x) || ischar(x));
ip.addParameter('behaviors', ...
    ["time", "vel", "lindist", "rewarded", "traj", "Y", ...
     "outBoundChoiceTimes", "inBoundChoiceTimes", "outBoundChoiceTimes", ...
     "rewardTimes", "errorTimes", "wellTimes", ...
     "pastRewarded", "pastLeftRight", "futureRewarded", "futureLeftRight", ...
     ],...
    @(x) iscellstr(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;
if ~isempty(Opt.figAppend) && ~endsWith(Opt.figAppend,"_")
    Opt.figAppend = Opt.figAppend + "_";
end
figAppend = Opt.figAppend + Option.animal;

% Determine columns of behavior that are logicals
logical_cols = arrayfun(@(x) islogical(behavior.(x)), Opt.behaviors);
logical_vars = Opt.behaviors(logical_cols);
% Determine columns of behavior with discrete, but not continuous values
discrete_cols = arrayfun(@(x) isnumeric(behavior.(x)) && ...
                    length(unique(behavior.(x))) < 10, Opt.behaviors);
discrete_vars = Opt.behaviors(discrete_cols);


%% Get the behavior
if isempty(behavior)
    running_times = r.timeBinMidPoints(r.sessionTypePerBin == 1);
    [behavior, throwout_times] = table.behavior.lookup(Option.animal, ...
        running_times);
end

%% Get the components

if isempty(Opt.names)
    if isfield(Option, "patternNamesFull")
        Opt.names = Option.patternNamesFull;
    else
        Opt.names = ["theta", "delta", "ripple"];
    end
end

if strcmpi(Opt.use, 'raw')
    activities = Components.activities;
elseif strcmpi(Opt.use, 'smooth')
    activities  = Components.smooth_activities;
else
    error("Invalid option for 'use': " + Opt.use);
end
time           = Components.time;


%% Plots





