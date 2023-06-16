% This script

%%
plots.plotOptDimHist;

%%

for p = 1:Option.numPartition
    for i = 1:nPatterns
        if floor(median_to_use(i)) == Patterns(p,1,i).rankRegress.optDimReducedRankRegress && 
            sameSource{i} = [sameSource{i}, Patterns(p,1,i)];
        end
        if floor(median_to_use(i)) == Patterns(p,2,i).rankRegress.optDimReducedRankRegress
            diffSource{i} = [diffSource{i}, Patterns(p,2,i)];
        end
    end
end

perf_same = cell(1,3);
perf_diff = cell(1,3);
for i = 1:nPatterns
    perf_same{i} = zeros(numel(sameSource{i}), floor(median_to_use(i)));
    perf_diff{i} = zeros(numel(diffSource{i}), floor(median_to_use(i)));
end

for i = 1:nPatterns
    num_same = size(perf_same{i},1);
    num_diff = size(perf_diff{i},1);
    for j = 1:num_same
        [perf_same{i}(j,:), ~] = plots.sequentialRemovePredDims(sameSource{i}(j).X_source, ...
            sameSource{i}(j).X_target, sameSource{i}(j).rankRegress.B_, sameSource{i}(j).rankRegress.optDimReducedRankRegress,...
            sameSource{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
    end
    for j = 1:num_diff
        [perf_diff{i}(j,:), ~] = plots.sequentialRemovePredDims(diffSource{i}(j).X_source, ...
            diffSource{i}(j).X_target, diffSource{i}(j).rankRegress.B_, diffSource{i}(j).rankRegress.optDimReducedRankRegress,...
            diffSource{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
    end
end