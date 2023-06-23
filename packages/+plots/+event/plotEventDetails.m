function plotEventDetails(Events, Option, varargin)
% plotEventDetails(Events, Option, varargin)
%
% Plots the details of the events in the Events structure.
%
% INPUTS
%   Events: structure containing the events
%   Option: structure containing the options
%   varargin: optional arguments
%       'saveFigures': whether to save the figures (default: false)
%       'appendFigTitle': string to append to the figure title
%       'savePath': path to save the figures
%       'lowerQuantile': lower quantile to use for the window
%       'upperQuantile': upper quantile to use for the window
%       'logAxis': which axes to log
%       'r': structure containing the spikes
%       'spikePlotStyle': 'heatmap' or 'raster' or 'heatmapBehind', 'rasterBehind'
%           heatmap : heatmap of spike counts
%           raster : raster plot of spikes
%           heatmapBehind : heatmap of spike counts behind the event time series
%           rasterBehind : raster plot of spikes behind the event time series
%       'showMuaSC': whether to show the MUA spike counts

disp("Plotting event details...");
tic

[Opt, H, nPatterns, times, windowQ, windowQlow] = ...
    plots.event.processEventsPreamble(Events, Option, varargin{:});
    % processEventsPreamble

%% EVENT DISTRIBUTIONS

% Histogram of magnitude and where we're setting high/low thresholds
f=figure();
set(f, 'Position', get(0, 'Screensize'));
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
    + newline + "WindowQuantile = " + Option.quantileToMakeWindows);
file = fullfile(Opt.savePath, ...
    "Histogram_of_event_magnitudes" + Opt.appendFigTitle + "_WindowQuantile=" + Option.quantileToMakeWindows)
savefig(f, file + ".fig");
saveas(f, file + ".png");
saveas(f, file + ".svg");

%% TIME SERIES

% Plot the events in the time series
f=fig("Time series of events" + Opt.appendFigTitle);
set(f, 'Position', get(0, 'Screensize'));
spikesGiven = ~isempty(Opt.r);
if spikesGiven && strcmp(Opt.spikePlotStyle, 'heatmap')
    tile = tiledlayout(nPatterns + 2, 1, ...
        'TileSpacing', 'compact', 'Padding', 'compact');
elseif spikesGiven && strcmp(Opt.spikePlotStyle, 'raster')
else
    tile = tiledlayout(nPatterns, 1, ...
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
        plots.heatmapBinSpikes(Opt.r, 'colorCols', 1:3, 'ylim', ylim());
    end
    if Opt.showMuaSC
        hold on
        plots.plotMuaSC(Opt.r, 'ylim', ylim());
    end
    ylabel(Option.patternNames(i) + " magnitude" + append);
end
linkaxes(tile.Children, 'x');

if ~isempty(Opt.r) && strcmp(Opt.spikePlotStyle, 'heatmap')

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
elseif ~isempty(Opt.r) && strcmp(Opt.spikePlotStyle, 'raster')
    error("Not implemented");
    % ax = nexttile();
    % plot(Opt.r);
    % title("Spike raster");
    % xlabel("Time (s)");
    % ylabel("Cell");
end
sgtitle("Time series of events" + Opt.appendFigTitle ...
    + newline + "WindowQuantile = " + Option.quantileToMakeWindows);
file = "Time_series_of_events" + Opt.appendFigTitle + "_WindowQuantile=" + Option.quantileToMakeWindows
file = fullfile(Opt.savePath, file);
savefig(file + ".fig");
saveas(gcf, file + ".png");
saveas(gcf, file + ".svg");


