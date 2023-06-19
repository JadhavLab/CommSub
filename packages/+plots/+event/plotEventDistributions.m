function plotEventDistributions(Events, Option, varargin)
% plotEventDistributions(Events, Option)
%
% Plots the distributions of the events in the Events structure.

disp("Plotting event distributions...");
tic


[Opt, H, nPatterns, times, windowQ, windowQlow] = ...
    plots.event.processEventsPreamble(Events, Option, varargin)

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
end

