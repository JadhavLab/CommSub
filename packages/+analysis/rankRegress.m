function Patterns = rankRegress(Patterns, Option)

if Option.waysOfPartitions ~= 2
    nTarget = size(Patterns(1).X_target,1);
    nSource = min(size(Patterns(1).X_source,1),...
        size(Patterns(1,2,1).X_source,1));
else
    nTarget = unique(arrayfun(@(x) numel(x.index_target), Patterns));
    nSource = unique(arrayfun(@(x) numel(x.index_source), Patterns));
end

% Number of cross validation folds.
numDimsUsedForPrediction = 1:min(nTarget, nSource);
cvNumFolds = Option.rankRegress.cvnum;
cvOptions = statset('crossval');
regressMethod = @ReducedRankRegress;
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) ...
    RegressFitAndPredict(regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
    numDimsUsedForPrediction, 'LossMeasure', 'NSE','RidgeInit', ...
    false, 'Scale', false);

% ------------------------------------------------------------------
% NOTES
% ------------------------------------------------------------------
% RegressMethod is what we use to apply
% RegressFitAndPredict has two lines: regression function and loss eval
% RankRegressRoutine : LIterally applies the crossval function to
% RegressFitAndPredict, and gives the cvLoss matrix and optimal dimension.
% and then recomputes B,V,B_ with optimal dimension
% -------------------------------------------------------------------

for n = progress(numel(Patterns), 'Title', 'RankRegress')

    p = Patterns(n);

    % when the partition is three-ways, j==1 means same target/source
    % pair and j==2 means diff target/source pair
    % disp("processing rrr for "+p+" partition and the "+i+" pattern 
    % "+j+" direction")
    curr_source = p.X_source';
    curr_target = p.X_target';

    nan_rows = any(isnan(curr_source), 2) | ...
               any(isnan(curr_target), 2); % setect

    [p.rankRegress.cvl, ...
     p.rankRegress.cvLoss, ...
     p.rankRegress.optDimReducedRankRegress,...
     p.rankRegress.B,...
     p.rankRegress.B_,...
     p.rankRegress.V] ...
        = rankRegressRoutine(cvFun, cvNumFolds, ...
              cvOptions, ...
              curr_target(~nan_rows, :), ...
              curr_source(~nan_rows,:), ...
              numDimsUsedForPrediction);

     p.rankRegress.muOptLoss   = p.rankRegress.cvLoss(1, p.rankRegress.optDimReducedRankRegress);
     p.rankRegress.stdOptLoss  = p.rankRegress.cvLoss(2, p.rankRegress.optDimReducedRankRegress);
     p.rankRegress.muFullLoss  = p.rankRegress.cvLoss(1, end);
     p.rankRegress.stdFullLoss = p.rankRegress.cvLoss(2, end);
    
    % Single neuron prediction
    if Option.analysis.singleNeuronPrediction

        B_singleprediction   = cell(1,nSource);
        dim_singleprediction = cell(1,nSource);

        for k = 1:nSource
            curr_singlesource = curr_source(:,k);
            if clean.zeroFiring(curr_singlesource)
                continue;
            end
            [~,~, dim_singleprediction{k}, B_singleprediction{k},~,~] = ...
                rankRegressRoutine(cvFun,...
                                   cvNumFolds, ...
                                   cvOptions, ...
                                   curr_target, ...
                                   curr_singlesource,...
                                   numDimsUsedForPrediction);
        end

        p.rankRegress.singlesource_B = B_singleprediction;
        p.rankRegress.singlesource_optDim = dim_singleprediction;
        p.rankRegress.B_rrr = getReducedB_(p.rankRegress.B,...
                                                     p.rankRegress.V, ...
                                                     nSource, ...
                                                     nTarget, ...
                                                     p.rankRegress.optDimReducedRankRegress);
    end


    Patterns(n) = p;
end
