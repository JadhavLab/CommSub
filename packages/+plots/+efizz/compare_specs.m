clear fields
clear fieldnames
% Define the frequency bands
freq_bands = struct('theta', [6, 10], 'delta', [0.5, 4], 'ripple', [150, 200]);

% Initialize structures to hold the mean time series and correlation coefficients
mean_time_series = struct('Cavg', struct('theta', [], 'delta', [], 'ripple', []), ...
                          'wpli', struct('theta', [], 'delta', [], 'ripple', []));
corr_coeffs      = struct('theta', [], 'delta', [], 'ripple', []);

% Iterate over each frequency band
F = fieldnames(freq_bands);
for i = 1:length(F)
    F = fieldnames{i};
    
    % Get the indices of the frequencies within the current frequency band
    freq_indices = find(efizz.f >= freq_bands.(freq_band)(1) & efizz.f <= freq_bands.(freq_band)(2));
    
    % Extract the time series for the current frequency band
    Cavg_time_series = mean(efizz.Cavg(:, freq_indices), 2);
    wpli_time_series = mean(efizz.wpli(:, freq_indices), 2);
    
    % Store the mean time series
    mean_time_series.Cavg.(freq_band) = Cavg_time_series;
    mean_time_series.wpli.(freq_band) = wpli_time_series;
    
    % Compute the correlation coefficient between the Cavg and wpli time series
    corr_coeffs.(freq_band) = corrcoef(Cavg_time_series, wpli_time_series);
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
