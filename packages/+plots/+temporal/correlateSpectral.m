function Components = correlateSpectal(Components, Events, Option, varargin)
% correlateSpectral(Components, Events, Option, varargin)
%   Correlate spectral components with events
%   
% Inputs:
%   Components - struct containing spectral components
%   Events     - struct containing events
%   Option     - struct containing options
%   varargin   - additional arguments
%
% Outputs:
%   none

for i = 1:numel(Components)
    Components(i).spec = singleCorrelateSpectral(Components(i), Events, Option, varargin{:});
end

% -------------------------------------------------------------------------

function out = singleCorrelateSpectral(Components, Events, Option, varargin)
% correlateSpectral(Components, Events, Option, varargin)
%   Correlate spectral components with events
% 
%  same as above, but for one single struct
%   

ip = inputParser();
ip.addParameter('names', [], @(x) iscellstr(x)  || isstring(x));
ip.addParameter('use', 'smooth', @(x) ischar(x) || isstring(x)); % 'smooth' or 'raw'
ip.addParameter('ci', 68, @(x) isscalar(x) && x > 0 && x < 100);
ip.addParameter('samples', 500, @(x) isscalar(x) && x > 0);
ip.addParameter('figAppend', "", @(x) ischar(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;

figAppend = Opt.animal + " " + Opt.figAppend;

%% Get components

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
time             = Components.time;


%% Plots

Htimes = Events.times;
Hvals  = Events.Hvals;
interpActivities = interp1(time, activities', Htimes);

% -------------------------------------------------------------------------

% PLOT: spearman correlation matrix between events and activities
combined = [interpActivities, Hvals];
corrMatrix = corr(combined, 'type', 'Spearman');
fig('spearman correlation matrix between events and activities')
imagesc(corrMatrix)
clim = [-1,1];
caxis(clim)
xticklabels(["comp1", "comp2", "comp3", Opt.names])
yticklabels(["comp1", "comp2", "comp3", Opt.names])
colormap(cmocean('balance'))
colorbar
out.spectral_correlation = corrMatrix;
savefig(gcf, fullfile(figuredefine(),"spectral correlation matrix" + figAppend + ".fig"))
saveas(gcf, fullfile(figuredefine(), "spectral correlation matrix" + figAppend + ".png"))
saveas(gcf, fullfile(figuredefine(), "spectral correlation matrix" + figAppend + ".svg"))

% -------------------------------------------------------------------------

% % PLOT: Cross-correlation between each event and activity
% % for each event, compute the cross-covariance between the event and each
% % activity component
% combined = [interpActivities, Hvals];
% labels = ["comp1", "comp2", "comp3", Opt.names];
% deltaT = median(diff(Htimes));
% fig('Cross-correlation between each event and activity components');clf
% ta = tiledlayout(size(combined,2), size(combined,2), 'TileSpacing', 'compact', 'Padding', 'compact');
% for i = 1:size(combined,2)
%     for j = 1:size(combined,2)
%         [labelA, labelB] = deal(labels(i), labels(j));
%         [xcov, lags] = xcorr(combined(:,i), combined(:,j), 100, 'coeff');
%         lagtimes = lags * deltaT;
%         ax = nexttile;
%         if i >= j
%             ax.Visible = 'off';
%         else
%             plot(lagtimes, xcov);
%             title("Cross-covariance between activity " + newline + labelA + " and " + labelB)
%         end
%     en_getripalignspiking_ver6tmp4.m:figure;plot(xcorr(mean(allSigVarPSTHsPFC(:,50:950)),mean(allSigVarPSTHsCA1(:,50:950))))
% cj_trodes2ff/Code/Src_Matlab/sj_HPExpt/HP_alluser/HPexpt_getfractionactiveSWRs.md
% end
% sgtitle("Cross-covariance between activity components")
% set(gcf, 'Position', [100, 100, 1000, 1000])

% -------------------------------------------------------------------------

% PLOT: CROSS-CORRELATION between each event and activity, subsampled
N = Opt.samples;
sampleSize = floor(size(combined,1)/N);
fig('Cross-correlation samples between each event and activity components');clf
ta = tiledlayout(10, 10, 'TileSpacing', 'compact', 'Padding', 'compact'); % Adjust layout according to your needs
% total = size(combined,2) * size(combined,2);
% cnt=0;
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
ci = Opt.ci;
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
out.cross_corr          = results;
out.cross_corr_mean     = means;
out.cross_corr_ci_lower = ci_lower;
out.cross_corr_ci_upper = ci_upper;
saveas(gcf, fullfile(figuredefine(), "cross_corr_" + figAppend + ".fig"))
saveas(gcf, fullfile(figuredefine(), "cross_corr_" + figAppend + ".png"))
saveas(gcf, fullfile(figuredefine(), "cross_corr_" + figAppend + ".svg"))

% -------------------------------------------------------------------------

% PLOT: Cross-covariance between each event and activity, subsampled
N = Opt.samples;
sampleSize = floor(size(combined,1)/N);
fig('Cross-covariance samples between each event and activity components');clf
ta = tiledlayout(10, 10, 'TileSpacing', 'compact', 'Padding', 'compact'); % Adjust layout according to your needs
% total = size(combined,2) * size(combined,2);
% cnt=0;
results = cell(size(combined,2), size(combined,2));
for k = progress(1:N, 'title', 'Subsampling')
    sampleStart = 1 + (k-1) * sampleSize;
    sampleEnd = k * sampleSize;
    subsample = combined(sampleStart:sampleEnd, :);
    for i = 1:size(subsample,2)
        for j = 1:size(subsample,2)
            [labelA, labelB] = deal(labels(i), labels(j));
            [xcov, lags] = xcov(subsample(:,i), subsample(:,j), min(300, sampleSize-1), 'coeff');
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
ci = Opt.ci;
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
out.cross_xcov          = results;
out.cross_xcov_mean     = means;
out.cross_xcov_ci_lower = ci_lower;
out.cross_xcov_ci_upper = ci_upper;
saveas(gcf, fullfile(figuredefine(), "cross_xcov_" + figAppend + ".fig"))
saveas(gcf, fullfile(figuredefine(), "cross_xcov_" + figAppend + ".png"))
saveas(gcf, fullfile(figuredefine(), "cross_xcov_" + figAppend + ".svg"))
