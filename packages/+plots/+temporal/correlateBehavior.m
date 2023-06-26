function correlateBehavior(Components, Behavior, varargin)
% correlateBehavior(Components, Behavior, varargin)
%
%   Correlate the behavior with the components of the model.
%
%   Inputs:

ip = inputParser();
ip.addParameter('behavior', []);
ip.addParameter('unique_times', []);
ip.addParameter('throwout_times', []);
ip.parse(varargin{:});
Opt = ip.Results;
figAppend = Opt.animal + " " + Opt.figAppend;

%% Get the behavior
if isempty(Opt.behavior)
    running_times = r.timeBinMidPoints(r.sessionTypePerBin == 1);
    [behavior, throwout_times] = table.behavior.lookup(Opt.animal, ...
        running_times);
    [behavior, unique_times] = behaviors.addBehToTable(behavior);
else
    behavior = Opt.behavior;
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
time             = Components.time;


%% Plots

