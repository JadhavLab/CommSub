function output = correlateBehavior(Components, Option, behavior, varargin)
% correlateBehavior(Components, Behavior, varargin)
%
%   Correlate the behavior with the components of the model.
%
%   Inputs:

ip = inputParser();
ip.addParameter('behavior', []);
ip.addParameter('unique_times', []);
ip.addParameter('throwout_times', []);
ip.addParameter('use', 'smooth');
ip.addParameter('names', []);
ip.addParameter('figAppend', "", @(x) isstring(x) || ischar(x));
ip.addParameter('behaviors', ...
    ["time", "vel", "lindist", "rewarded", "traj", "Y", ...
     "outBoundChoiceTimes", "inBoundChoiceTimes", "outBoundChoiceTimes", ...
     "rewardTimes", "errorTimes", "wellTimes", ...
     "pastRewarded", "pastLeftRight", "futureRewarded", "futureLeftRight", ...
     ],...
    @(x) iscellstr(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;
if ~isempty(Opt.figAppend) && ~endsWith(Opt.figAppend,"_")
    Opt.figAppend = Opt.figAppend + "_";
end
figAppend = Opt.figAppend + Option.animal;

% Determine columns of behavior that are logicals
logical_cols = arrayfun(@(x) islogical(behavior.(x)), Opt.behaviors);
logical_vars = Opt.behaviors(logical_cols);
% Determine columns of behavior with discrete, but not continuous values
discrete_cols = arrayfun(@(x) isnumeric(behavior.(x)) && ...
                    length(unique(behavior.(x))) < 10, Opt.behaviors);
discrete_vars = Opt.behaviors(discrete_cols);

folder="corrBehavior";

%% Get the behavior
if isempty(behavior)
    running_times = r.timeBinMidPoints(r.sessionTypePerBin == 1);
    [behavior, throwout_times] = table.behavior.lookup(Option.animal, ...
        running_times);
end

%% Get the components
if isempty(Opt.names)
    if isfield(Option, "patternNamesFull")
        Opt.names = Option.patternNamesFull;
    else
        Opt.names = ["theta", "delta", "ripple"];
    end
end

if strcmpi(Opt.use, 'raw')
    activities = Components.activities;
elseif strcmpi(Opt.use, 'smooth')
    activities  = Components.smooth_activities;
else
    error("Invalid option for 'use': " + Opt.use);
end
time           = Components.time;

% Interp
Btimes = behavior.time;
behavior = behavior(:, Opt.behaviors);
Bvals    = table2array(behavior(:, Opt.behaviors));
interpActivities = interp1(time, activities', Btimes);

%% Plots

% PLOT: spearman correlation matrix between behaviors and activities
combined = [interpActivities, Bvals];
corrMatrix = corr(combined, 'type', 'Spearman');
fig('spearman correlation matrix between behaviors and activities')
imagesc(corrMatrix)
clim = [-1,1];
caxis(clim)
xticklabels([Opt.names, "Behavior"])
yticklabels([Opt.names, "Behavior"])
colormap(cmocean('balance'))
colorbar
savefig(gcf, fullfile(figuredefine(folder),"behavior correlation matrix" + figAppend + ".fig"))
saveas(gcf, fullfile(figuredefine(folder), "behavior correlation matrix" + figAppend + ".png"))
saveas(gcf, fullfile(figuredefine(folder), "behavior correlation matrix" + figAppend + ".svg"))


%% PLOT: 1 Cross-correlation between each behavior and activity
% for each behavior, compute the cross-covariance between the behavior and each
% activity component
combined = [interpActivities, Bvals];
labels = [Opt.names, "Behavior"];
deltaT = median(diff(Btimes));
fig('Cross-correlation between each behavior and activity components');clf
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
sgtitle("Cross-covariance between activity components and behaviors")
set(gcf, 'Position', [100, 100, 1000, 1000])


%% PLOT: 2 Time course of the behavior, the activity, and the cross-covariance
fig('Time course of the behavior, the activity, and the cross-covariance')
for i = 1:size(activities, 2)
    subplot(size(activities, 2), 1, i)
    plot(time, activities(:, i), 'b');
    hold on
    plot(Btimes, Bvals, 'r');
    legend([Opt.names(i), 'Behavior'])
    [xcov, lags] = xcorr(activities(:, i), Bvals, 'coeff');
    plot(time, xcov, 'k');
    legend('Activity', 'Behavior', 'Cross-covariance')
    xlabel('Time')
    ylabel('Activity/Behavior')
    title(['Time course of activity ', Opt.names(i), ', behavior, and their cross-covariance'])
    grid on
end
savefig(gcf, fullfile(figuredefine(folder),"time course behavior activity" + figAppend + ".fig"))
saveas(gcf, fullfile(figuredefine(folder), "time course behavior activity" + figAppend + ".png"))
saveas(gcf, fullfile(figuredefine(folder), "time course behavior activity" + figAppend + ".svg"))

%% PLOT: 3 Moving window cross-covariance
nSubsamples = 100; % replace with your value
subsampleSize = floor(length(time) / nSubsamples); % replace with your value

% Initialize a 3D matrix to hold the correlation results
Rvals = nan(nSubsamples, length(behaviorVars), length(componentNames));

% Initialize a 3D matrix to hold the p-value results
Pvals = nan(nSubsamples, length(behaviorVars), length(componentNames));

% Loop over subsamples
for k = 1:nSubsamples
    % Get the indices for this subsample
    idxStart = (k - 1) * subsampleSize + 1;
    idxEnd = min(k * subsampleSize, length(time)); % Use min in case length(time) is not a multiple of subsampleSize

    % Extract the subsample of activities and time
    subsampleActivities = activities(:, idxStart:idxEnd);
    subsampleTime = time(idxStart:idxEnd);

    % Loop over behavior variables and components
    for i = 1:length(behaviorVars)
        for j = 1:length(componentNames)
            % Extract the subsample of behavior
            subsampleBehavior = behaviors.(behaviorVars{i})(idxStart:idxEnd);

            % Calculate the correlation and p-value for this behavior and component
            [Rvals(k, i, j), Pvals(k, i, j)] = corr(subsampleActivities(j, :)', subsampleBehavior, 'type', 'Spearman');
        end
    end
end



%% PLOT: 4 Permutation test
numPermutations = 1000;
permutationXcovs = zeros(numPermutations, size(activities, 2));
for i = 1:size(activities, 2)
    [actualXcov, ~] = xcorr(activities(:, i), Bvals, 'coeff');
    for j = 1:numPermutations
        shuffledBvals = Bvals(randperm(length(Bvals)));
        permutationXcovs(j, i) = xcorr(activities(:, i), shuffledBvals, 'coeff');
    end
    permpValue(i) = mean(permutationXcovs(:, i) > actualXcov);
    disp(['p-value for activity ', Opt.names(i), ': ', num2str(pValue)])
end

%% PLOT: 5 Granger causality
%[F, pVal] =  plots.temporal.grangercausality(activities, Bvals, maxlag);
grang = struct();
for i = 1:size(Bvals, 2)
    [grang.F, grang.pVal, grang.issues] = ... 
        plots.temporal.grangerCausality(activities(:, 1:3), Bvals(:, i), 100);
    disp(['p-value for activity ', Opt.names(i), ': ', num2str(pVal)]);
    % Compare to last 3 components
    [grang.last_F, grang.last_pVal, grang.last_issues] = ...
        plots.temporal.grangerCausality(activities(:, end-2:end), Bvals(:, i), 100);
    disp(['p-value for private activity ', Opt.names(i), ': ', num2str(last_pVal)]);
    % Compare to shuffled behavior
    [grang.shuffled_F, grang.shuffled_pVal, grant.shuff_issues] = ...
        plots.temporal.grangerCausality(activities(:, 1:3), Bvals(randperm(length(Bvals)), i), 100);
end

% Collect the p-values and F statistics for each analysis
pVals = [grang.pVal; grang.last_pVal; grang.shuffled_pVal]';
Fs = [grang.F; grang.last_F; grang.shuffled_F]';

% Create a figure
figure

% Create a subplot for p-values
subplot(2,1,1)
bar(pVals)
set(gca, 'XTickLabel',Opt.names)
legend('First 3 components','Last 3 components','Shuffled')
ylabel('p-value')
title('Granger Causality p-values')
set(gca,'YScale','log') % Optional: use a logarithmic scale to make small p-values easier to compare

% Create a subplot for F statistics
subplot(2,1,2)
bar(Fs)
set(gca, 'XTickLabel',Opt.names)
legend('First 3 components','Last 3 components','Shuffled')
ylabel('F statistic')
title('Granger Causality F statistics')

% Store computed output values from entire script
output.grang = grang;
output.permutationXcovs = permutationXcovs;
output.permutationPvals = permpValue;
output.xcorr.Rvals = Rvals;
output.xcorr.Pvals = Pvals;
output.xcorr.lags = lags;
output.corrMatrix = corrMatrix;
