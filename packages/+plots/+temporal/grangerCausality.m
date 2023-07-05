function [F, pVal, issues] = grangerCausality(activities, Bvals, maxLag)
%GRANGERCAUSALITY Perform Granger causality test on the given data
%   [F, pVal] = GRANGERCAUSALITY(activities, Bvals, maxLag) performs a
%   Granger causality test on the given data. The test is performed by
%   regressing the B-values on their own lagged values and the lagged
%   values of the activities. The F-statistic and p-value of the test are
%   returned.
%
%   Inputs:
%   - activities: a vector of activity values
%   - Bvals: a vector of B-values
%   - maxLag: the maximum lag to use in the test
%
if ~adftest(activities) || ~adftest(Bvals)
    issues.stationarity = false;
else
    issues.stationarity = true;
end

% Create matrix of lagged values
activities_lagged = lagmatrix(activities, 1:maxLag);
Bvals_lagged = lagmatrix(Bvals, 1:maxLag);

% Remove rows with NaN values (created by lagging process)
validIdx = all(~isnan(activities_lagged), 2) & all(~isnan(Bvals_lagged), 2);
activities_lagged = activities_lagged(validIdx, :);
Bvals_lagged = Bvals_lagged(validIdx, :);

% Perform regression of Bvals on its own lagged values and the lagged values of activities
X = [ones(size(Bvals_lagged, 1), 1), Bvals_lagged, activities_lagged];
mdl_full = fitlm(X, Bvals(validIdx));

% Check for autocorrelation in residuals
if lbqtest(mdl_full.Residuals.Raw)
    issues.autocorr_mdl_full = true;
else
    issues.autocorr_mdl_full = false;
end

% Perform regression of Bvals on its own lagged values only
X = [ones(size(Bvals_lagged, 1), 1), Bvals_lagged];
mdl_reduced = fitlm(X, Bvals(validIdx));

% Check for autocorrelation in residuals
if lbqtest(mdl_reduced.Residuals.Raw)
    issues.autocorr_mdl_reduced = true;
else
    issues.autocorr_mdl_reduced = false;
end

% Calculate the F-statistic
F = ((mdl_reduced.SSE - mdl_full.SSE) / mdl_full.NumEstimatedCoefficients) / mdl_full.MSE;

% Determine the p-value of the F-statistic
pVal = 1 - fcdf(F, mdl_full.NumEstimatedCoefficients, mdl_full.DFE);
