function [] = runChecks(Events, Spk, Patterns, Option)
% Run checks on the data overall for raw
% and processed data
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

cellOfWindows = Events.cellOfWindows;

disp("Running checks")
% firing rate checks
a = @() plots.frChecks(Spk, "appendFigTitle", char(Option.animal));
b = @() plots.event.plotEventDetails(Events, Option, ...
    "savePath", char(fullfile(figuredefine(), "eventDetails")), ...
    "appendFigTitle", char(Option.animal), ...
    'Spk', Spk, ...
    'spikePlotStyle', 'heatmapBehind');
c = @() plots.event.plotSingleEvents(cellOfWindows, Events, Option, Spk, ...
    'before', 0.5, 'after', 0.5, 'appendFigTitle', char(Option.animal), ...
    'eventTypes', ["theta", "delta", "ripple"], ...
    'displayType', 'heatmap');
d = @() plots.event.plotSingleEvents(cellOfWindows, Events, Option, Spk, ...
    'before', 0.5, 'after', 0.5, 'appendFigTitle', char(Option.animal), ...
    'eventTypes', ["theta", "delta", "ripple"], ...
    'displayType', 'raster');
% figure
% Run plots in parallel
a = parfeval(gcp, a, 0);
b = parfeval(gcp, b, 0);
c = parfeval(gcp, c, 0);
d = parfeval(gcp, d, 0);
Wait for plots to finish
wait(a); wait(b); wait(c);
% a(); b(); c();

