function Patterns = JPECC(Patterns, Option)
% Computes lagged JPECC between patterns in two brain araes

% Number of cross validation folds.
cvNumFolds = Option.jpecc.cvnum;

cvOptions = statset('crossval');
regressMethod = @ReducedRankRegress;
numDimsUsedForPrediction = 1:min(nTarget,nSource);
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
    (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
    numDimsUsedForPrediction, 'LossMeasure', 'NSE','RidgeInit', ...
    false, 'Scale', false);

nTarget = size(P(1).X_target,1);
nSource = min(size(P(1).X_source,1), size(P(1,2,1).X_source,1));

% ------------------------------------------------------------------
% NOTES
% ------------------------------------------------------------------
% RegressMethod is what we use to apply
% RegressFitAndPredict has two lines: regression function and loss eval
% RankRegressRoutine : LIterally applies the crossval function to
% RegressFitAndPredict, and gives the cvLoss matrix and optimal dimension.
% and then recomputes B,V,B_ with optimal dimension
% -------------------------------------------------------------------

% Loop and apply
for i = progress(1:numel(P), 'Title', 'Rank regression')

    p = P(i);
        
    curr_source = p.X_source(:, :)';
    curr_target = p.X_target(:, :)';
    tens_curr_source = p.X_source;
    tens_curr_target = p.X_target;

    nan_rows = any(isnan(curr_source), 2) | ...
               any(isnan(curr_target), 2);

    [p.rankRegress.cvl,...
     p.rankRegress.cvLoss,...
     p.rankRegress.optDimReducedRankRegress,...
     p.rankRegress.B, ...
     p.rankRegress.B_, ...
     p.rankRegress.V] ...
        = rankRegressRoutine(cvFun, cvNumFolds, ...
          cvOptions, curr_target(~nan_rows, :),...
          curr_source(~nan_rows, :), ...
          numDimsUsedForPrediction);

     p.rankRegress.muOptLoss   = p.rankRegress.cvLoss(1, p.rankRegress.optDimReducedRankRegress);
     p.rankRegress.stdOptLoss  = p.rankRegress.cvLoss(2, p.rankRegress.optDimReducedRankRegress);
     p.rankRegress.muFullLoss  = p.rankRegress.cvLoss(1, end);
     p.rankRegress.stdFullLoss = p.rankRegress.cvLoss(2, end);

    % JPECC : Check all time lags
    for iTime = 1:size(tens_curr_source, 2)
    for jTime = 1:size(tens_curr_source, 2)

        cs = squeeze(tens_curr_source(:, iTime, :))';
        ct = squeeze(tens_curr_target(:, jTime, :))';
        nan_rows = any(isnan(cs), 2) | any(isnan(ct), 2);

        lessThan3Samples = sum(~nan_rows) < 3;
        if lessThan3Samples
            p.jpecc(iTime, jTime).cvl    = nan;
            p.jpecc(iTime, jTime).cvLoss = nan;
            p.jpecc(iTime, jTime).optDimReducedRankRegress = nan;
            p.jpecc(iTime, jTime).muOptLoss   = nan;
            p.jpecc(iTime, jTime).stdOptLoss  = nan;
            p.jpecc(iTime, jTime).muFullLoss  = nan;
            p.jpecc(iTime, jTime).stdFullLoss = nan;
            continue
        end

        [p.jpecc(iTime, jTime).cvl,...
         p.jpecc(iTime, jTime).cvLoss,...
         p.jpecc(iTime, jTime).optDimReducedRankRegress] ...
            = rankRegressRoutine(cvFun, cvFolds, cvOptions, ct(~nan_rows, :), cs(~nan_rows, :), numDimsUsedForPrediction);
         p.jpecc(iTime, jTime).muOptLoss   = p.jpecc(iTime, jTime).cvLoss(1, p.rankRegress.optDimReducedRankRegress);
         p.jpecc(iTime, jTime).stdOptLoss  = p.jpecc(iTime, jTime).cvLoss(2, p.rankRegress.optDimReducedRankRegress);
         p.jpecc(iTime, jTime).muFullLoss  = p.jpecc(iTime, jTime).cvLoss(1, end);
         p.jpecc(iTime, jTime).stdFullLoss = p.jpecc(iTime, jTime).cvLoss(2, end);

    end
    end

    % Single neuron prediction
    if doSingleNeuronPrediction
        % not implemented for jpecc
    end

    P(i) = p;

end
