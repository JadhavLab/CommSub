
out = struct();

source = Patterns.X_source;
target = Patterns.X_target;
time   = Patterns.X_time;
source_index = Patterns.index_source;
target_index = Patterns.index_target;
source = r.spikeCountMatrix(source_index,:);
target = r.spikeCountMatrix(target_index,:);
time = r.timeBinMidPoints;

rr = Patterns.rankRegress
B_ = rr.B_;
[u,s,v] = svd(B_, 'full');
sdiag = diag(s);

disp("Size source: " + size(source))
disp("Size u: " + size(u))
disp("Size target: " + size(target))
disp("Size v: " + size(v))

% ----------------
% OLD METHOD :: can FP when activity in one area but not the other
% ----------------
% s_u = s(1:3,1:3);
% s_v = s(1:3,1:3);
s_u = 1;
s_v = 1;
both = [source;target];
match_inputoutput = [u(:,1:3)*s_u;v(:,1:3)*s_v];
size(match_inputoutput), size(both)
oldactivities = match_inputoutput' * (both);
oldactivities = activities - mean(activities,2);
oldactivities = abs(activities);
out.oldactivities = oldactivities;


% ----------------
% NEW METHOD :: can't FP when activity in one area but not the other
% ----------------
activation_source = u(:,1:3)'*source;
activation_target = v(:,1:3)'*target;
activities = activation_source .* activation_target;
smooth_activties = zeros(size(activities));
for i = progress(1:size(activities,1))
    smooth_activties(i,:) = smooth(activities(i,:), 400);
end
out.activities = activities;
out.smooth_activties = smooth_activties;

% We're doing the above now in analysis.temporal.match()

figure
subplot(2,1,1)
plot(abs(activities'))
legend({'component 1', 'component 2', 'component 3'})
title('Matched input-output components')
subplot(2,1,2)
plot(abs(smooth_activties'))
legend({'component 1', 'component 2', 'component 3'})
title('Smoothed matched input-output components')

% PLOT: single events
plots.event.plotSingleEvents(cellOfWindows, Events, Option, r, ...
    'before', 0.5, 'after', 0.5, 'appendFigTitle', char(Option.animal), ...
    'eventTypes', ["theta", "delta", "ripple"], ...
    'subspaceComponents', activities, ...
    'subspaceTimes', time, ...
    'saveFolder', fullfile(figuredefine, 'singleEvent_plusActivities'), ...
    'displayType', 'heatmap');

% PLOT: correlation with events
Htimes = Events.times;
Hvals = Events.Hvals;
% interpolate activity to event times
interpActivities = interp1(time, smooth_activties', Htimes);

% spearman correlation matrix between events and activities
combined = [interpActivities, Hvals];
corrMatrix = corr(combined, 'type', 'Spearman');
figure
imagesc(corrMatrix)
clim = [-1,1];
caxis(clim)
xticklabels(["comp1", "comp2", "comp3", "theta", "delta", "ripple"])
yticklabels(["comp1", "comp2", "comp3", "theta", "delta", "ripple"])
colormap(cmocean('balance'))
colorbar

% PLOT: Cross-covariance between each event and activity
% for each event, compute the cross-covariance between the event and each
% activity component
combined = [interpActivities, Hvals];
labels = ["comp1", "comp2", "comp3", "theta", "delta", "ripple"];
deltaT = median(diff(Htimes));
figure
ta = tiledlayout(size(combined,2), size(combined,2), 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:size(combined,2)
    for j = 1:size(combined,2)
        [labelA, labelB] = deal(labels(i), labels(j));
        [xcov, lags] = xcorr(combined(:,i), combined(:,j), 100, 'coeff');
        lagtimes = lags * deltaT;
        ax = nexttile;
        if i >= j
            ax.Visible = 'off';
        else
            plot(lagtimes, xcov);
            title("Cross-covariance between activity " + newline + labelA + " and " + labelB)
        end
    end
end
sgtitle("Cross-covariance between activity components")
set(gcf, 'Position', [100, 100, 1000, 1000])

% PLOT: Cross-covariance between each event and activity, subsampled
N = 500;
sampleSize = floor(size(combined,1)/N);
figure
ta = tiledlayout(10, 10, 'TileSpacing', 'compact', 'Padding', 'compact'); % Adjust layout according to your needs
total = size(combined,2) * size(combined,2);
cnt=0;
results = cell(size(combined,2), size(combined,2));
for k = progress(1:N, 'title', 'Subsampling')
    sampleStart = 1 + (k-1) * sampleSize;
    sampleEnd = k * sampleSize;
    subsample = combined(sampleStart:sampleEnd, :);
    for i = 1:size(subsample,2)
        for j = 1:size(subsample,2)
            [labelA, labelB] = deal(labels(i), labels(j));
            [xcov, lags] = xcorr(subsample(:,i), subsample(:,j), min(300, sampleSize-1), 'coeff');
            lagtimes = lags * deltaT;
            results{i,j} = [results{i,j}, xcov];
            cnt = sub2ind([size(combined,2), size(combined,2)], i, j);
            ax = subplot(size(combined,2), size(combined,2), cnt);
            if i >= j
                ax.Visible = 'off';
            else
                hold on
                plot(ax, lagtimes, xcov, 'Color', [0.5,0.5,0.5,0.5], ...
                                     'LineWidth', 0.5, 'LineStyle', ':');
                alpha(0.33)
                if k == 1
                    title(ax, labelA + newline + labelB)
                    ylim([-1,1])
                    if j == 1
                        ylabel(ax, labelA + newline + "Cross-covariance")
                    end
                    if i == size(combined,2)
                        xlabel(ax, labelB + newline + "Lag (s)")
                    end
                end
            end
        end
    end
end
sgtitle("Cross-covariance between activity components across subsamples")
set(gcf, 'Position', [100, 100, 1000, 1000])
% determine the mean, and ci% confidence interval
disp("Computing mean and confidence interval")
ci = 0.68;
means = cellfun(@(x) mean(x,2), results, 'UniformOutput', false);
ci_lower = cellfun(@(x) prctile(x, ci/2, 2), results, 'UniformOutput', false);
ci_upper = cellfun(@(x) prctile(x, 1-(ci/2), 2), results, 'UniformOutput', false);
for i = progress(1:size(combined,2), 'title', 'Plotting means and confidence intervals')
    for j = 1:size(combined,2)
        ax = subplot(size(combined,2), size(combined,2), ...
             sub2ind([size(combined,2), size(combined,2)], i, j));
        if i < j
            hold on
            plot(ax, lagtimes, means{i,j}, 'Color', 'k', 'LineWidth', 2);
            plot(ax, lagtimes, ci_lower{i,j}, 'Color', 'b', 'LineWidth', 1);
            plot(ax, lagtimes, ci_upper{i,j}, 'Color', 'b', 'LineWidth', 1);
        end
    end
end

