clear fields
clear fieldnames
% Define the frequency bands
freq_bands = struct('theta', [6, 10], 'delta', [0.5, 4], 'ripple', [150, 200]);

% Initialize structures to hold the mean time series and correlation coefficients
mean_time_series = struct('Cavg', struct('theta', [], 'delta', [], 'ripple', []), ...
                          'wpli', struct('theta', [], 'delta', [], 'ripple', []),...
                            'S1', struct('theta', [], 'delta', [], 'ripple', []), ...
                            'S2', struct('theta', [], 'delta', [], 'ripple', []));
                    
corr_coeffs      = struct('theta', [], 'delta', [], 'ripple', []);

% Iterate over each frequency band
F = fieldnames(freq_bands);
for i = 1:length(F)
    
    f = F{i};
    % Get the indices of the frequencies within the current frequency band
    freq_indices = find(efizz.f >= freq_bands.(f)(1) & efizz.f <= freq_bands.(f)(2));
    % Extract the time series for the current frequency band
    Cavg_time_series = mean(efizz.Cavg(:, freq_indices), 2);
    wpli_time_series = mean(efizz.wpli(:, freq_indices), 2);
    s1_time_series   = mean(efizz.S1(:, freq_indices), 2);
    s2_time_series   = mean(efizz.S2(:, freq_indices), 2);
    % Store the mean time series
    mean_time_series.Cavg.(f) = Cavg_time_series;
    mean_time_series.wpli.(f) = wpli_time_series;
    mean_time_series.S1.(f)   = s1_time_series;
    mean_time_series.S2.(f)   = s2_time_series;
    % Compute the correlation coefficient between the Cavg and wpli time series
    corr_coeffs.(f) = corrcoef(Cavg_time_series, wpli_time_series);
end

% % And cross-correlation coefficients
% for i = 1:length(F)
% for j = 1:length(F)
%     corr_coeffs.(F{i} + "_" + F{j}) = corrcoef(mean_time_series.Cavg.(F{i}), mean_time_series.wpli.(F{j}));
% end
% end

% Plot the time series
figure;
subplot(2, 1, 1);
hold on;
plot(mean_time_series.Cavg.theta,  'r');
plot(mean_time_series.Cavg.delta,  'g');
plot(mean_time_series.Cavg.ripple, 'b');
title('Mean Time Series for Cavg');
legend('Theta', 'Delta', 'Ripple');
hold off;

subplot(2, 1, 2);
hold on;
plot(mean_time_series.wpli.theta,  'r');
plot(mean_time_series.wpli.delta,  'g');
plot(mean_time_series.wpli.ripple, 'b');
title('Mean Time Series for wpli');
legend('Theta', 'Delta', 'Ripple');
hold off;

linkaxes(get(gcf, 'children'), 'x');

% Print the correlation coefficients
disp('Correlation Coefficients:');
disp(struct2table(corr_coeffs));

% Correlation of ripple coherence/wpli to ripple S1 and ripple S2 band
C = corrcoef([mean_time_series.wpli.ripple, mean_time_series.S1.ripple, mean_time_series.S2.ripple]);

% Load global ripple time
ndb.load(Option.animal, "globalripple");
dat = [];
globalripple = globalripple{1};
for i = 1:length(globalripple)
    if isempty(globalripple{i})
        continue
    end
    dat = [dat; globalripple{i}];
end
% dat = [start, stop , amp] matrix

% get mean_time_series.wpli.ripple inside versus outside of start stop ripple
% windows
wpli = mean_time_series.wpli.ripple;
t = efizz.t;
in = []; % mean wpli inside ripple windows
out = []; % mean wpli outside ripple windows
tstart = interp1(t, 1:length(t), dat(:,1), 'nearest');
tstop  = interp1(t, 1:length(t), dat(:,2), 'nearest');
out = [out, mean(wpli(1:tstart(1)-1))];
for i = 1:length(dat(:,1))
    in = [in, mean(wpli(tstart(i):tstop(i+1)))];
    if i > 1
        out = [out, mean(wpli(tstop(i-1)+1:tstart(1)-1))];
    end
end
out = [out, mean(wpli(tstop(end)+1:end))];

figure
% bootstrap_in_mean = bootstrp(1000, @nanmean, in);
% bootstrap_out_mean = bootstrp(1000, @nanmean, out);
ci_in = bootci(1000, @nanmean, in);
ci_out = bootci(1000, @nanmean, out);
bar(1:2, [nanmean(in), nanmean(out)]); hold on;
errorbar(1:2, [nanmean(in), nanmean(out)], [nanmean(in)-ci_in(1), nanmean(out)-ci_out(1)], [ci_in(2)-nanmean(in), ci_out(2)-nanmean(out)], 'k.');
xticks([1, 2]);
xticklabels({'Inside', 'Outside'});
title('Mean wpli inside vs outside ripple windows');
