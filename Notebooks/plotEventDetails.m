function plotEventDetails(Events, Option, varargin)
% plotEventDetails(Events)
%
% Plots the details of the events in the Events structure.

disp("Plotting event details...");
tic

% Parse input
ip = inputParser();
ip.addParameter('saveFigures', false, @islogical);
ip.addParameter('appendFigTitle', '', @ischar);
ip.addParameter('savePath', '', @ischar);
ip.addParameter('lowerQuantile', 0.001, @isnumeric);
ip.addParameter('upperQuantile', 0.99, @isnumeric);
ip.addParameter('logAxis', [2], @isnumeric); % log delta axis
ip.addParameter('r', [], @isstruct); % r structure -- if given, we plot the 
                                     % spikes
ip.addParameter('spikePlotStyle', 'heatmap', @ischar); 
% 'heatmap' or 'raster' or 'heatmapBehind', 'rasterBehind'
% heatmap : heatmap of spike counts
% raster : raster plot of spikes
% heatmapBehind : heatmap of spike counts behind the event time series
% rasterBehind : raster plot of spikes behind the event time series
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

%% EVENT DISTRIBUTIONS

% Histogram of magnitude and where we're setting high/low thresholds
figure();
tile = tiledlayout(nPatterns, 1, ...
    'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:nPatterns
    nexttile();
    histogram(H(:,i), 'NumBins', 1000, 'Normalization', 'count')
    title(Option.patternNames(i));
    hold on;
    line([windowQ(i) windowQ(i)], ylim(), 'Color', 'r');
    line([windowQlow(i) windowQlow(i)], ylim(), 'Color', 'r', 'LineStyle', ':');
    % Plot the CDF
    [N,edges] = histcounts(H(:,i), 'Normalization', 'cdf', 'NumBins', 1000 );
    bin_centers = edges(1:end-1) + diff(edges)/2;
    N = N * max(ylim());
    plot(bin_centers, N, 'b:', 'LineWidth', 1.5); % light black dashed line
    if i == nPatterns
        legend("", "Window quantile (high)", "Window quantile (low)");
    end
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
f=fig("Time series of events" + Opt.appendFigTitle);
figure(f)
spikesGiven = ~isempty(Opt.r)
if spikesGiven && strcmp(Opt.spikePlotStyle, 'heatmap')
    tile = tiledlayout(nPatterns, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
elseif spikesGiven && strcmp(Opt.spikePlotStyle, 'raster')
else
    tile = tiledlayout(nPatterns + 2, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
end
for i = 1:nPatterns
    nexttile();
    h = H(:,i);
    if ismember(i, Opt.logAxis); h = log(h); end
    plot(times, h);
    title(Option.patternNames(i));
    hold on;
    upper = [windowQ(i) windowQ(i)];
    lower = [windowQlow(i) windowQlow(i)];
    if ismember(i, Opt.logAxis)
        upper = log(upper);
        lower = log(lower);
    end
    line([times(1) times(end)], upper, 'Color', 'r');
    line([times(1) times(end)], lower, 'Color', 'r', 'LineStyle', ':');
    if i == nPatterns
        legend("", "Window quantile (high)", "Window quantile (low)");
        xlabel("Time (s)");
    end
    append = "";
    if ismember(i, Opt.logAxis)
        append = " (log)";
    end
    if Opt.spikePlotStyle == "heatmapBehind"
        hold on
        plots.heatmapBinSpikes(r, 'colorCols', 1:3, 'ylim', ylim());
    end
    ylabel(Option.patternNames(i) + " magnitude" + append);
end
linkaxes(tile.Children, 'x');
if ~isempty(Opt.r)

    scm = Opt.r.spikeCountMatrix;
    hpc = Opt.r.celllookup.region == "CA1";
    pfc = Opt.r.celllookup.region == "PFC";
    scm_hpc = scm(hpc, :);
    scm_pfc = scm(pfc, :);
    scm_hpc(scm_hpc == 0) = nan;
    scm_pfc(scm_pfc == 0) = nan;
    scm_times = Opt.r.timeBinMidPoints;

    % ca1
    ax = nexttile();
    cmap = gray;
    % cmap = flipud(cmap);
    cmap(:, 1:2) = 0;
    cmap = cmap(128:end, :);
    cmap = repelem(cmap, 2, 1);
    cmap(1, :) = 1;
    % cmap(1, :) = 1;
    imagesc(scm_times, 1:size(scm_hpc, 1), scm_hpc);
    colormap(ax,cmap)
    colorbar();
    title("Spike rate");
    xlabel("Time (s)");
    ylabel("Spike rate");

    % pfc
    ax = nexttile();
    cmap = gray;
    % cmap = flipud(cmap);
    cmap(:, 2:3) = 0;
    cmap = cmap(128:end, :);
    cmap = repelem(cmap, 2, 1);
    cmap(1, :) = 1;
    imagesc(scm_times, 1:size(scm_pfc, 1), scm_pfc);
    colormap(ax,cmap)
    colorbar();
end

