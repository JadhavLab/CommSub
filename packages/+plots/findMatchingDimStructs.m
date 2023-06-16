plots.plotOptDimHist;

%% go on and run the three-way version

waysOfPartitions = 3;
Option.targetArea = "PFC";
 
clear Patterns_ReSplitThreeWays
Patterns_ReSplitThreeWays = struct("X_source",[], "X_target",[]);
Patterns_ReSplitThreeWays.rankRegress = struct(...
    "B", [], ...
    "B_", [], ...
    "optDimReducedRankRegress", 0);

patternNames = ["theta","delta","ripple",...
    "theta-control","delta-control","ripple-control"];


Patterns_ReSplitThreeWays = repmat(Patterns_ReSplitThreeWays, [Option.numPartition,2,numResult]);
for iPartition = 1:Option.numPartition
    % check which is the target region
    if Option.targetArea == "CA1"
        [same_target, diff_target, t, st_index, dt_index, t_index ] ...
            = trialSpikes.threeGroupsSplit(hpcFR, pfcFR, X_hpc, X_pfc, Option.binsToMatchFR);
    else
        [same_target, diff_target, t, st_index, dt_index, t_index ] ...
            = trialSpikes.threeGroupsSplit(pfcFR, hpcFR, X_pfc, X_hpc, Option.binsToMatchFR);
    end
    for i = 1:nPatterns*2
        Patterns_ReSplitThreeWays(iPartition,1,i).X_source = same_target{i};
        Patterns_ReSplitThreeWays(iPartition,2,i).X_source = diff_target{i};
        Patterns_ReSplitThreeWays(iPartition,1,i).X_target = t{i};
        Patterns_ReSplitThreeWays(iPartition,2,i).X_target = t{i};
        
        Patterns_ReSplitThreeWays(iPartition,1,i).index_source = st_index;
        Patterns_ReSplitThreeWays(iPartition,2,i).index_source = dt_index;
        Patterns_ReSplitThreeWays(iPartition,1,i).index_target = t_index;
        Patterns_ReSplitThreeWays(iPartition,2,i).index_target = t_index;
        
        if Option.targetArea == "CA1"
            Patterns_ReSplitThreeWays(iPartition,1,i).directionality = "hpc-hpc";
            Patterns_ReSplitThreeWays(iPartition,2,i).directionality = "pfc-hpc";
        else
            Patterns_ReSplitThreeWays(iPartition,1,i).directionality = "pfc-pfc";
            Patterns_ReSplitThreeWays(iPartition,2,i).directionality = "hpc-pfc";
        end
        
        Patterns_ReSplitThreeWays(iPartition,1,i).name = patternNames(i);
        Patterns_ReSplitThreeWays(iPartition,2,i).name = patternNames(i);
    end
end


%%
%%%%%%%%%%%%%%%% RANK-REGRESS SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if waysOfPartitions ~= 2
    nTarget = size(Patterns_ReSplitThreeWays(1,1,1).X_target,1);
    nSource = min(size(Patterns_ReSplitThreeWays(1,1,1).X_source,1),...
        size(Patterns_ReSplitThreeWays(1,2,1).X_source,1));
end
numDimsUsedForPrediction = 1:min(nTarget,nSource);

for p = 1:Option.numPartition
    for i = 1:numResult
        
        B_singleprediction = cell(1,nSource);
        dim_singleprediction = cell(1,nSource);
        
        % Number of cross validation folds.
        cvNumFolds = 10;
        cvOptions = statset('crossval');
        regressMethod = @ReducedRankRegress;
        cvFun = @(Ytrain, Xtrain, Ytest, Xtest) RegressFitAndPredict...
            (regressMethod, Ytrain, Xtrain, Ytest, Xtest, ...
            numDimsUsedForPrediction, 'LossMeasure', 'NSE','RidgeInit', ...
            false, 'Scale', false);
        
        
        % when the partition is three-ways, j==1 means same target/source
        % pair and j==2 means diff target/source pair
        for j = [HPC, PFC]
            curr_source = (Patterns_ReSplitThreeWays(p,j,i).X_source)';
            curr_target = (Patterns_ReSplitThreeWays(p,j,i).X_target)';
            [   Patterns_ReSplitThreeWays(p,j,i).rankRegress.cvl, ...
                Patterns_ReSplitThreeWays(p,j,i).rankRegress.cvLoss,...
                Patterns_ReSplitThreeWays(p,j,i).rankRegress.optDimReducedRankRegress,...
                Patterns_ReSplitThreeWays(p,j,i).rankRegress.B,...
                Patterns_ReSplitThreeWays(p,j,i).rankRegress.B_,...
                Patterns_ReSplitThreeWays(p,j,i).rankRegress.V] ...
                = rankRegressRoutine(cvFun, cvNumFolds, ...
                cvOptions, curr_target, curr_source, ...
                numDimsUsedForPrediction);
            Patterns_ReSplitThreeWays(p,j,i).rankRegress.B_rrr = getReducedB_(Patterns_ReSplitThreeWays(p,j,i).rankRegress.B,...
                Patterns_ReSplitThreeWays(p,j,i).rankRegress.V, nSource, nTarget,...
                Patterns_ReSplitThreeWays(p,j,i).rankRegress.optDimReducedRankRegress);
        end
    end
end

%% find the matching structs with the median

same_source = Patterns_ReSplitThreeWays(:,1,:);
diff_source = Patterns_ReSplitThreeWays(:,2,:);

same = [];
same_theta = [];
same_delta = [];
same_ripple = [];

diff = [];
diff_theta = [];
diff_delta = [];
diff_ripple = [];

for p = 1:Option.numPartition
    for i = 1:nPatterns
        if same_source(p,1,i).rankRegress.optDimReducedRankRegress == floor(median_for_each_pattern(i))
            same = [same, same_source(p,1,i)];
            if (i==1)
                same_theta = [same_theta, numel(same)];
            elseif (i==2)
                same_delta = [same_delta, numel(same)];
            else
                same_ripple = [same_ripple, numel(same)];
            end
        end
        if diff_source(p,1,i).rankRegress.optDimReducedRankRegress == floor(median_for_each_pattern(i))
            diff = [diff, diff_source(p,1,i)];
            if (i==1)
                diff_theta = [diff_theta, numel(diff)];
            elseif (i==2)
                diff_delta = [diff_delta, numel(diff)];
            else
                diff_ripple = [diff_ripple, numel(diff)];
            end
        end
    end
end

%% remove for same/diff per pattern type

perf_same = [];
perf_same_theta = [];
perf_same_delta = [];
perf_same_ripple = [];


perf_diff = [];
perf_diff_theta = [];
perf_diff_delta = [];
perf_diff_ripple = [];

for i = 1:numel(same)
    
          [same(i).rankRegress.removedPerformance, ~,~]=...
                        plots.getUncorrelatedPerformance( same(i).rankRegress.B_,...
                        same(i).X_source, same(i).X_target, same(i).rankRegress.optDimReducedRankRegress,...
                        numDimsUsedForPrediction, same(i).rankRegress.cvLoss);
    
           [same(i).perfAsRemovingEach] = plots.sequentialRemovePredDims(same(i).X_source, ...
                    same(i).X_target,same(i).rankRegress.B_, same(i).rankRegress.optDimReducedRankRegress, ...
                    same(i).rankRegress.cvLoss,numDimsUsedForPrediction);
    if (ismember(i, same_theta))
        perf_same_theta = [perf_same_theta; same(i).perfAsRemovingEach];
    elseif (ismember(i, same_delta))
        perf_same_delta = [perf_same_delta; same(i).perfAsRemovingEach];
    else 
        perf_same_ripple = [perf_same_ripple; same(i).perfAsRemovingEach];
    end
%     perf_same = [perf_same;same(i).perfAsRemovingEach];
end
%%

for i = 1:numel(diff)
          [diff(i).rankRegress.removedPerformance, ~,~]=...
                        plots.getUncorrelatedPerformance( diff(i).rankRegress.B_,...
                        diff(i).X_source, diff(i).X_target, diff(i).rankRegress.optDimReducedRankRegress,...
                        numDimsUsedForPrediction, diff(i).rankRegress.cvLoss);
    
           [diff(i).perfAsRemovingEach] = plots.sequentialRemovePredDims(diff(i).X_source, ...
                    diff(i).X_target,diff(i).rankRegress.B_, diff(i).rankRegress.optDimReducedRankRegress, ...
                    diff(i).rankRegress.cvLoss,numDimsUsedForPrediction);
                
    if (ismember(i, diff_theta))
        perf_diff_theta = [perf_diff_theta; diff(i).perfAsRemovingEach];
    elseif (ismember(i, diff_delta))
        perf_diff_delta = [perf_diff_delta; diff(i).perfAsRemovingEach];
    else 
        perf_diff_ripple = [perf_diff_ripple; diff(i).perfAsRemovingEach];
    end
%     perf_diff = [perf_diff; diff(i).perfAsRemovingEach];
end
%%
% calcualte std and make error bars out of them
% 
% perf_same_mean = mean(perf_same,1);
% perf_diff_mean = mean(perf_diff, 1);
% perf_same_std = std(perf_same,1);
% perf_diff_std = std(perf_diff,1);


perf_same_mean_theta = mean(perf_same_theta,1);
perf_diff_mean_theta = mean(perf_diff_theta, 1);
perf_same_std_theta = std(perf_same_theta,1);
perf_diff_std_theta = std(perf_diff_theta,1);


perf_same_mean_delta = mean(perf_same_delta,1);
perf_diff_mean_delta = mean(perf_diff_delta, 1);
perf_same_std_delta = std(perf_same_delta,1);
perf_diff_std_delta = std(perf_diff_delta,1);


perf_same_mean_ripple = mean(perf_same_ripple,1);
perf_diff_mean_ripple = mean(perf_diff_ripple, 1);
perf_same_std_ripple = std(perf_same_ripple,1);
perf_diff_std_ripple = std(perf_diff_ripple,1);
%%
figure(667)


% % plot(perf_same_mean);
% errorbar(perf_same_mean, perf_same_std)
% hold on
% % plot(perf_diff_mean)
% errorbar(perf_diff_mean, perf_diff_std);
% 
% legend("predicting same region", "predicting different region")
% xlabel("number of preditive dimensions removed")
% ylabel("performance")

subplot(1,3,1)
% plot(perf_same_mean);
if numel(perf_same_std_theta) <=1
    plot(perf_same_mean_theta);
else
errorbar(perf_same_mean_theta, perf_same_std_theta);
end
% errorbar(perf_same_mean_theta, perf_same_std_theta)
hold on
% plot(perf_diff_mean)
if numel(perf_diff_std_theta) <=1
   plot(perf_diff_mean_theta);
else
errorbar(perf_diff_mean_theta, perf_diff_std_theta);
end

legend("predicting same region", "predicting different region")
xlabel("number of preditive dimensions removed")
ylabel("performance theta")

subplot(1,3,2) 
% plot(perf_same_mean);
if numel(perf_same_std_delta) <=1
    plot(perf_same_mean_delta);
else
errorbar(perf_same_mean_delta, perf_same_std_delta);
end

% errorbar(perf_same_mean_delta, perf_same_std_delta)
hold on
% plot(perf_diff_mean)
if numel(perf_diff_std_delta) <=1
    plot(perf_diff_mean_delta);
else
errorbar(perf_diff_mean_ delta, perf_diff_std_delta);
end

legend("predicting same region", "predicting different region")
xlabel("number of preditive dimensions removed")
ylabel("performance delta")

subplot(1,3,3)
% plot(perf_same_mean);
if numel(perf_same_std_ripple) <=1
    plot(perf_same_mean_ripple);
else
errorbar(perf_same_mean_ripple, perf_same_std_ripple);
end
% errorbar(perf_same_mean_ripple, perf_same_std_ripple)
hold on
% plot(perf_diff_mean)
if numel(perf_diff_std_ripple) <=1
   plot(perf_diff_mean_ripple);
else
errorbar(perf_diff_mean_ripple, perf_diff_std_ripple);
end

legend("predicting same region", "predicting different region")
xlabel("number of preditive dimensions removed")
ylabel("performance ripple")

