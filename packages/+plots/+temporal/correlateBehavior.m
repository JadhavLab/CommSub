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
ip.addParameter('use', 'smooth');
ip.addParameter('names', []);
ip.addParameter('figAppend', "", @(x) isstring(x) || ischar(x));
ip.parse(varargin{:});
Opt = ip.Results;
if ~isemtpy(Opt.figAppend) && ~endsWith(Opt.figAppend,"_")
    Opt.figAppend = Opt.figAppend + "_";
end
end
figAppend = Opt.figAppend + Option.animal;

%% Get the behavior
if isempty(Opt.behavior)
    running_times = r.timeBinMidPoints(r.sessionTypePerBin == 1);
    [behavior, throwout_times] = table.behavior.lookup(Option.animal, ...
        running_times);
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
time           = Components.time;


%% Plots

