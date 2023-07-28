function out = regressefizz(efizz, Patterns_overall, field, varargin)
    % Input parser
    p = inputParser;
    addParameter(p, 'faxis', inf, @isnumeric);
    addParameter(p, 'tabPrepend', "");
    addParameter(p, 'ploton', false, @islogical);
    parse(p, varargin{:});
    Opt = p.Results;
    
    % Get the downsampling factor
    d_factor = p.Results.faxis;

    % Extract time series
    efizz_t = efizz.t;
    Patterns_overall_t = Patterns_overall.X_time;

    % If downsampling is required in frequency domain
    if ~isempty(d_factor) && d_factor > 1 && d_factor ~= inf
        % Get the number of samples in the original data
        original_samples = size(efizz.(field), 2);

        % Create a vector of original sample points
        original_sample_points = linspace(min(efizz.f), max(efizz.f), original_samples);
        % Create a vector of desired sample points
        desired_sample_points = linspace(min(efizz.f), max(efizz.f), d_factor);
        % Interpolate the data to the desired number of sample points
        if strcmp(field, 'phi')
            % If the field is 'phi', we interpolate real and imaginary parts separately
            real_part = interp1(original_sample_points, real(exp(1i*efizz.(field))), desired_sample_points);
            imag_part = interp1(original_sample_points, imag(exp(1i*efizz.(field))), desired_sample_points);
            efizz.(field) = angle(real_part + 1i*imag_part);
        else
            % For other fields, we can interpolate directly
            tmp = zeros(size(efizz.(field), 1), d_factor);
            for i = progress(1:size(efizz.(field), 1))
                tmp(i, :) = interp1(original_sample_points, efizz.(field)(i, :), desired_sample_points);
            end
            efizz.(field) = tmp;
        end
    else
        desired_sample_points = efizz.f;
    end

    % Interpolate to align the time series
    aligned_efizz = interp1(efizz_t, efizz.(field), Patterns_overall_t, 'nearest', 'extrap');

    % Extract the first 5 canonical vectors
    U = Patterns_overall.cca.u(:, 1:5);
    V = Patterns_overall.cca.v(:, 1:5);

    % dropnan 
    nancols = any(isnan(aligned_efizz), 1);
    aligned_efizz(:, nancols) = 0;
    nanlocs = any(isnan(aligned_efizz), 2) | any(isnan(U), 2) | any(isnan(V), 2);
    aligned_efizz(nanlocs, :) = [];
    U(nanlocs, :) = [];
    V(nanlocs, :) = [];

    for i = 1:5

        if field ~= "phi"
            % Perform regression for U
            model_U = fitlm(aligned_efizz, U(:, i));
            % Perform regression for V
            model_V = fitlm(aligned_efizz, V(:, i));

            % Get coefficients
            coef_U = model_U.Coefficients.Estimate;
            coef_V = model_V.Coefficients.Estimate;

            % Get F-statistics and p-values
            F_U = model_U.Coefficients.tStat;
            pvalue_U = model_U.Coefficients.pValue;
            F_V = model_V.Coefficients.tStat;
            pvalue_V = model_V.Coefficients.pValue;
        else
            % Convert to unit circle coordinates
            x = cosd(aligned_efizz);
            y = sind(aligned_efizz);
            % Perform regression for U
            model_U_x = fitlm(U, x);
            model_U_y = fitlm(U, y);
            % Perform regression for V
            model_V_x = fitlm(V, x);
            model_V_y = fitlm(V, y);
            % Get coefficients
            coef_U_x = model_U_x.Coefficients.Estimate;
            coef_U_y = model_U_y.Coefficients.Estimate;
            coef_V_x = model_V_x.Coefficients.Estimate;
            coef_V_y = model_V_y.Coefficients.Estimate;
            coef_U = sqrt(coef_U_x^2 + coef_U_y^2);
            coef_V = sqrt(coef_V_x^2 + coef_V_y^2);
            % Get F-statistics and p-values
            F_U_x = model_U_x.Coefficients.tStat;
            F_U_y = model_U_y.Coefficients.tStat;
            F_V_x = model_V_x.Coefficients.tStat;
            F_V_y = model_V_y.Coefficients.tStat;
            F_U = sqrt(F_U_x^2 + F_U_y^2);
            F_V = sqrt(F_V_x^2 + F_V_y^2);
            pvalue_U_x = model_U_x.Coefficients.pValue;
            pvalue_U_y = model_U_y.Coefficients.pValue;
            pvalue_V_x = model_V_x.Coefficients.pValue;
            pvalue_V_y = model_V_y.Coefficients.pValue;
            % combine p-values with souffers method
            % pcomb = (1-erf(sum(sqrt(2) * erfinv(1-2*p))/sqrt(2*length(p))))/2;
            pvalue_U = (1-erf(sum(sqrt(2) * erfinv(1-2*[pvalue_U_x, pvalue_U_y]))/sqrt(2*2)))/2;
            pvalue_V = (1-erf(sum(sqrt(2) * erfinv(1-2*[pvalue_V_x, pvalue_V_y]))/sqrt(2*2)))/2;
        end

        out(i).F_U = F_U(2:end);
        out(i).pvalue_U = pvalue_U(2:end);
        out(i).F_V = F_V(2:end);
        out(i).pvalue_V = pvalue_V(2:end);
        out(i).coef_U = coef_U(2:end);
        out(i).coef_V = coef_V(2:end);
        out(i).coef_i = repmat(i, numel(F_U(2:end)), 1);
        out(i).f = desired_sample_points(:);
        out(i).field = repmat(string(field), numel(F_U(2:end)), 1);
    end

    tab = struct2table(out(1), 'AsArray', false);
    Opt.tabPrepend = Opt.tabPrepend + "_" + field;
    Opt.tabPrepend = Opt.tabPrepend + "_faxis=" + Opt.faxis;
    writetable(tab, figuredefine('tables', Opt.tabPrepend + "_regress.csv"));

    if Opt.ploton
        fig(Opt.tabPrepend + "_" + field + "_regress_faxis=" + Opt.faxis);
        plots.cca.plot_coefs(out, 0.01, 'figAppend', Opt.tabPrepend+"_regress");
    end


end

