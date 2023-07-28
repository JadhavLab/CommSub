function plot_coefs(out, efizz, alpha)
    if nargin < 3
        alpha = 0.05;
    end

    % Number of subplots
    num_plots = numel(out);

    % Extract frequencies
    freqs = efizz.f;

    for i = 1:num_plots
        % Create subplot
        subplot(num_plots, 1, i);

        % Get coefficients and p-values
        coef_U = out(i).coef_U;
        coef_V = out(i).coef_V;
        pvalue_U = out(i).pvalue_U;
        pvalue_V = out(i).pvalue_V;

        % Replace non-significant coefficients with NaN
        coef_U(pvalue_U > alpha) = NaN;
        coef_V(pvalue_V > alpha) = NaN;

        % Stem and scatter plot for coefficients
        stem(freqs, coef_U(2:end), 'b');
        hold on;
        % scatter(freqs, coef_U(2:end), 'b', 'filled', 'SizeData', 4);
        stem(freqs, coef_V(2:end), 'r');
        % scatter(freqs, coef_V(2:end), 'r', 'filled', 'SizeData', 4);
        hold off;

        % Set title and labels
        title(sprintf('Canonical Vector %d', i));
        xlabel('Frequency');
        ylabel('Coefficient');
        legend('U', 'V');

        sz = get(0, 'Screensize');
        % make 1/10 of width screen
        sz = [sz(3)+10, sz(4)+10, sz(3)/10, sz(4)];
        set(gcf, 'Position', sz);
    end
    sgtitle("Canonical Coefficients: " + out(1).name);
end

