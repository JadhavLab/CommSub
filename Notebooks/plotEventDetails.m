function plotEventDetails(Events, Option)
% plotEventDetails(Events)
%
% Plots the details of the events in the Events structure.

disp("Plotting event details...");
tic

ip = inputParser();
ip.addParameter('saveFigures', false, @islogical);
ip.addParameter('appendFigTitle', '', @ischar);
ip.addParameter('savePath', '', @ischar);
ip.addParameter('lowerQuantile', 0.001, @isnumeric);
ip.addParameter('upperQuantile', 0.99, @isnumeric);
ip.parse(varargin{:});
Opt = ip.Results;

if isempty(Opt.savePath)
    Opt.savePath = fullfile(figuredefine(), "EventDetails");
    if ~isfolder(Opt.savePath)
        mkdir(Opt.savePath);
    end
    disp("Save path for figures is " + Opt.savePath);
end


H = Events.H;
nPatterns = size(Events.H, 2);
quantileToMakeWindows = Option.quantileToMakeWindows;

% Nan above and below the quantiles
windowQ = zeros(nPatterns, 1);
for i = 1:nPatterns
    lowerQ = quantile(H(:, i), Opt.lowerQuantile);
    upperQ = quantile(H(:, i), Opt.upperQuantile);
    windowQ(i) = quantile(H(:, i), quantileToMakeWindows);
    disp("Pattern: " + Option.patternNames(i));
    disp("Lower quantile: " + lowerQ);
    disp("Upper quantile: " + upperQ);
    H(H(:, i) < lowerQ, i) = nan;
    H(H(:, i) > upperQ, i) = nan;
end
Hc = num2cell(H, 1);
assert(~isequal(Hc{:}), "All values are equal")

%% EVENT DISTRIBUTIONS

% Histogram of magnitude and where we're setting high/low thresholds
figure();
tile = tiledlayout(nPatterns, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nPatterns
    nexttile();
    histogram(H(:,i), 'BinWidth', 0.1);
    title(Option.patternNames(i));
    hold on;
    line([windowQ(i) windowQ(i)], ylim(), 'Color', 'r');
end
sgtitle("Histogram of event magnitudes" + Opt.appendFigTitle ...
    + newline + "WindowQuantile = " + quantileToMakeWindows);
savefig("Histogram_of_event_magnitudes" + Opt.appendFigTitle ...
    + "_WindowQuantile=" + quantileToMakeWindows + ".fig");
saveas(gcf, "Histogram_of_event_magnitudes" + Opt.appendFigTitle ...
    + "_WindowQuantile=" + quantileToMakeWindows + ".png");
saveas(gcf, "Histogram_of_event_magnitudes" + Opt.appendFigTitle ...
    + "_WindowQuantile=" + quantileToMakeWindows + ".svg");

%% TIME SERIES

% Plot the events in the time series

