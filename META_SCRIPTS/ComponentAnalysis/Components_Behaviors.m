% predict without any partition
% run stanard script up to the X matrices

nTarget = size(X_pfc{1},1);
nSource = size(X_hpc{1},1);

numDimsUsedForPrediction = 1:min(nTarget,nSource);
cvNumFolds = 10;
cvOptions = statset('crossval');
regressMethod = @ReducedRankRegress;
cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
    (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
    numDimsUsedForPrediction, 'LossMeasure', 'NSE','RidgeInit', ...
    false, 'Scale', false);
B_ = cell(1,6);
optDim = cell(1,6);
for i = 1:numel(patternNames)
    curr_target = X_pfc{i}';
    curr_source = X_hpc{i}';
    
    [~, ~,optDim{i},~,B_{i},~] ...
        = rankRegressRoutine(cvFun, cvNumFolds, cvOptions, curr_target,...
        curr_source, numDimsUsedForPrediction);
end

%% total subspaces for all patterns
%% --------------------------------
dim = 5; % chose the top five dimensions
celllookup = cellInfo.getCellIdentities(animal, cell_index, areaPerNeuron);
[total_subspaces, cell_subspaces] = components.subspaceSimilarity(dim, B_, spikeRateMatrix, celllookup,"CA1");


%% TIME-partitioned CALCULATIONS
%[animal_behavior, behavior_running_subspaces] = components.animalBehComponents...
%               (animal,timeBinMidPoints,sessionTypePerBin, total_subspaces);
[animal_behavior, throw_out_times] = table.behavior.lookup(animal, timeBinMidPoints, [], 'throwOutSleep', true);
running_subspaces = total_subspaces(:,~throwout_times);

%% make structs for the critical behaviors for plotting purposes
%components.makeBehaviorStructs
[critical_behaviors, unifiedTime, critical_components, critical_times] = ...
    components.makeBehaviorStructs(5, cellOfWindows, animal_behavior, false, 0.85);

