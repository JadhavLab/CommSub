function [] = runChecks(Option, Events, r, cellOfWindows)
% Run checks on the data overall for raw
% and processed data
%
% Parameters:
% -----------
% Option: struct
%   Options for the analysis
% Events: struct
%   Events for the analysis
% r: struct
%   Raw data
% cellOfWindows: cell
%   Windows for the analysis
%
% Returns:
% --------

disp("Running checks")
% firing rate checks
a = @() plots.frChecks(r, "appendFigTitle", char(Option.animal));
b = @() plots.event.plotEventDetails(Events, Option, ...
    "savePath", char(fullfile(figuredefine(), "eventDetails")), ...
    "appendFigTitle", char(Option.animal), ...
    'r', r, ...
    'spikePlotStyle', 'heatmapBehind');
c = @() plots.event.plotSingleEvents(cellOfWindows, Events, Option, r, ...
    'before', 0.5, 'after', 0.5, 'appendFigTitle', char(Option.animal), ...
    'eventTypes', ["theta", "delta", "ripple"], ...
    'displayType', 'heatmap');
% Run plots in parallel
% a = parfeval(gcp, a, 0);
% b = parfeval(gcp, b, 0);
% c = parfeval(gcp, c, 0);
% Wait for plots to finish
% wait(a); wait(b); wait(c);
a(); b(); c();

