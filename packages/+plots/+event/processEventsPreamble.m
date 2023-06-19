function [Opt, H, nPatterns, times, windowQ, windowQlow] = processEventsPreamble(Events, Option, varargin)
% processEventsPreamble(Events, Option)
%
% Handles the common preamble in processing event structures.

disp("Processing events preamble...");
tic

% Parse input
ip = inputParser();
ip.addParameter('saveFigures', false, @islogical);
ip.addParameter('appendFigTitle', '', @ischar);
ip.addParameter('savePath', '', @ischar);
ip.addParameter('lowerQuantile', 0.001, @isnumeric);
ip.addParameter('upperQuantile', 0.99, @isnumeric);
ip.addParameter('logAxis', [2], @isnumeric); % log delta axis
ip.addParameter('r', [], @isstruct); % r structure -- if given, we plot the spikes
ip.addParameter('spikePlotStyle', 'heatmapBehind', @ischar); 
ip.addParameter('showMuaSC', true, @islogical); % show MUA spike counts
ip.parse(varargin{:});
Opt = ip.Results;

if isempty(Opt.savePath)
    Opt.savePath = fullfile(figuredefine(), "EventDetails");
    if ~isfolder(Opt.savePath)
        mkdir(Opt.savePath);
    end
    disp("Save path for figures is " + Opt.savePath);
end

times = Events.times;
H = Events.Hvals;
nPatterns = size(Events.H, 2);
quantileToMakeWindows = Option.quantileToMakeWindows;
quantileToMakeWindows_low = 1 - quantileToMakeWindows;

% Nan above and below the quantiles
windowQ = zeros(nPatterns, 1);
windowQlow = zeros(nPatterns, 1);
for i = 1:nPatterns
    lowerQ = quantile(H(:, i), Opt.lowerQuantile);
    upperQ = quantile(H(:, i), Opt.upperQuantile);
    windowQ(i) = quantile(H(:, i), quantileToMakeWindows);
    windowQlow(i) = quantile(H(:, i), quantileToMakeWindows_low);
    disp("Pattern: " + Option.patternNames(i));
    disp("Lower quantile: " + lowerQ);
    disp("Upper quantile: " + upperQ);
    H(H(:, i) < lowerQ, i) = nan;
    H(H(:, i) > upperQ, i) = nan;
end
Hc = num2cell(H, 1);
assert(~isequal(Hc{:}), "All values are equal")

end

