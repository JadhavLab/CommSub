function plotTimeSeriesOfEvents(Events, Option, varargin)
% plotTimeSeriesOfEvents(Events, Option)
%
% Plots the time series of the events in the Events structure.

disp("Plotting time series of events...");
tic


[Opt, H, nPatterns, times, windowQ, windowQlow] = ...
    plots.event.processEventsPreamble(Events, Option, varargin)

%% TIME SERIES
% Plot the events in the time series
f=figure("Time series of events" + Opt.appendFigTitle);
figure(f)
spikesGiven = ~isempty(Opt.r)
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
        plots.heatmapBinSpikes(r, 'colorCols', 1:3, 'ylim', ylim());
    end
    if Opt.showMuaSC
        hold on
        plots.plotMuaSC(r, 'ylim', ylim());
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
    + newline + "WindowQuantile = " + quantileToMakeWindows);
savefig("Time_series_of_events" + Opt.appendFigTitle ...
    + "_WindowQuantile=" + quantileToMakeWindows + ".fig");

