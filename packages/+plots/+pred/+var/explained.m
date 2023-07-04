function [var_explained, mean_var_explained, nan_indices] = ...
                   calculateVarianceExplained(X_source, X_target, B)
% calculate the predictive performance of source firing to target firing
% outputs the performance and the mean performance or each neuron
% predicting the other 

[~,nTarget] = size(X_target);
try
    [~, yhat] = RegressPredict(X_target, X_source, B);
catch
    keyboard % error in prediction
end

    keyboard

var_explained = [];

for j = 1:nTarget
    unaccounted_variance = nansum((yhat(:,j)-X_target(:,j)).^2);
    total_variance       = nansum((nanmean(X_target(:,j))-X_target(:,j)).^2);
    temp_r_square = 1 - unaccounted_variance/total_variance;
    var_explained = [var_explained temp_r_square];
end

nan_indices = find(isnan(var_explained));
indices     = intersect(find(var_explained>0), find(~isnan(var_explained)));

% only take the non -Inf and non nan performances
mean_var_explained = nanmean(var_explained(indices));

end

