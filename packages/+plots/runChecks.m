function [] = runChecks(Events, Spk, Patterns, Option, varargin)
% Run checks on the data overall for raw and processed data
%
% Parameters:
% -----------
% Option: struct
%   Options for the analysis
% Events: struct
%   Events for the analysis
% Patterns: struct
%   Patterns for the analysis
% Spk: struct
%   Raw data
% cellOfWindows: cell
%   Windows for the analysis
%
% Returns:
% --------

ip = inputParser();
ip.addParameter('wait', true, @islogical);
ip.addParameter('parallel', true, @islogical);
ip.addParameter('visible', 'on', @ischar);
ip.parse(varargin{:});
Opt = ip.Results;

tic

cellOfWindows = Events.cellOfWindows;

disp("Running checks")
% firing rate checks
a = @() plots.frChecks(Spk, "appendFigTitle", char(Option.animal),...
    'visible', Opt.visible);
b = @() plots.event.plotEventDetails(Events, Option, ...
    "savePath", char(fullfile(figuredefine(), "eventDetails")), ...
    "appendFigTitle", char(Option.animal), ...
    'Spk', Spk, ...
    'spikePlotStyle', 'heatmapBehind', 'visible', Opt.visible);
c = @() plots.event.plotSingleEvents(cellOfWindows, Events, Option, Spk, ...
    'before', 0.5, 'after', 0.5, 'appendFigTitle', char(Option.animal), ...
    'eventTypes', ["theta", "delta", "ripple"], ...
    'displayType', 'heatmap', 'visible', Opt.visible);
% d = @() plots.event.plotSingleEvents(cellOfWindows, Events, Option, Spk, ...
%     'before', 0.5, 'after', 0.5, 'appendFigTitle', char(Option.animal), ...
%     'eventTypes', ["theta", "delta", "ripple"], ...
%     'displayType', 'raster');
% figure

if Opt.parallel
    disp("Parallel processing - enabled")
else
    disp("Parallel processing - disabled")
end

func_list = {a, b, c, d};
out_list = cell(1, length(func_list));
func_names = cellfun(@func2str, func_list, 'UniformOutput', false);

% Run plots in parallel
for i = 1:length(func_list)
    disp("Running " + func2str(func_list{i}));
    if Opt.parallel
        out_list{i} = parfeval(gcp, func_list{i}, 0);
    else
        out_list{i} = func_list{i}();
    end
end

% Wait for plots to finish
if Opt.parallel
    if Opt.wait
        disp("Waiting for plots to finish")
        for i = 1:length(func_list)
            disp("Waiting for " + func2str(func_list{i}));
            out_list{i}.wait();
        end
    end
    % Print any error messages
    disp("-------")
    disp("Errors")
    disp("-------")
    for i = 1:length(func_list)
        if ~isempty(out_list{i}.Error)
            disp(func_names{i} + " errors: " + out_list{i}.Error.message);
        end
    end
end

disp("Done running checks" + newline + toc + " seconds elapsed")

function out= d()
    out=struct();
    fig("event timeseries"); events.plot_events(Events, 'rewardTimes');
    fig("iei"); events.plotIEI(Events, 'NumBins', 1000);
    fig("cross iei"); events.plotCrossIEI(Events, 'NumBins', 1000, 'layout', 'matrix');
end


end % of runChecks

o
    
