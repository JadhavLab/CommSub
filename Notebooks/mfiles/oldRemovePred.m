% find median optimal dimensions for each pattern

[median_optDim, all_optdim, curr_pattern_dims, median_for_each_pattern, mean_for_each_pattern]...
    = plots.plotOptDimHist(Patterns, Option, numResult);
median_for_each_pattern = ceil(median_for_each_pattern);
mean_for_each_pattern   = ceil(mean_for_each_pattern); % direction x pattern
%% 
% 

% select optimal for all directions
median_to_use = zeros(1,nPatterns);
median_to_use = max(median_for_each_pattern, [], 1);

% select from the Patterns those fit the optdim
sameSource = cell(1,3);
diffSource = cell(1,3);

all_optdim = true; % whether to use a global median for

if all_optdim
    for p = 1:Option.numPartition
        for i = 1:nPatterns
            if floor(median_optDim) == Patterns(p,1,i).rankRegress.optDimReducedRankRegress && ...
                    floor(median_optDim) == Patterns(p,2,i).rankRegress.optDimReducedRankRegress
                sameSource{i} = [sameSource{i}, Patterns(p,1,i)];
                diffSource{i} = [diffSource{i}, Patterns(p,2,i)];
            end
        end
    end
else
    % not fully implemented yet - unsure how to make jagged arrays
    for p = 1:Option.numPartition
        for i = 1:nPatterns
            sameSource{i} = [sameSource{i}, Patterns(p,1,i)];
            diffSource{i} = [diffSource{i}, Patterns(p,2,i)];
        end
    end
end

%%
% make structures

perf_same_same = cell(1,3); % using the B from HPC-HPC to predict HPC
perf_same_diff = cell(1,3); % using the B from HPC-HPC to predict PFC

perf_diff_same = cell(1,3); % using the B from HPC-PFC to predict HPC
perf_diff_diff = cell(1,3); % using the B from HPC-PFC to predict PFC

for i = 1:nPatterns
    if all_optdim
        perf_same_same{i} = zeros(numel(sameSource{i}), floor(median_optDim));
        perf_same_diff{i} = zeros(numel(diffSource{i}), floor(median_optDim));
        perf_diff_diff{i} = zeros(numel(diffSource{i}), floor(median_optDim));
        perf_diff_same{i} = zeros(numel(diffSource{i}), floor(median_optDim));
    else
        perf_same_same{i} = zeros(numel(sameSource{i}), floor(median_to_use(i)));
        perf_same_diff{i} = zeros(numel(diffSource{i}), floor(median_to_use(i)));
        perf_diff_diff{i} = zeros(numel(diffSource{i}), floor(median_to_use(i)));
        perf_diff_same{i} = zeros(numel(diffSource{i}), floor(median_to_use(i)));
    end
    
end

for i = 1:nPatterns
    num_same = size(perf_same_same{i},1);
    for j = 1:num_same
        [perf_same_same{i}(j,:), ~] = plots.sequentialRemovePredDims(sameSource{i}(j).X_source, ...
            sameSource{i}(j).X_target, sameSource{i}(j).rankRegress.B_, sameSource{i}(j).rankRegress.optDimReducedRankRegress,...
            sameSource{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        [perf_same_diff{i}(j,:), ~] = plots.sequentialRemovePredDims(sameSource{i}(j).X_source, ...
            diffSource{i}(j).X_target, sameSource{i}(j).rankRegress.B_, sameSource{i}(j).rankRegress.optDimReducedRankRegress,...
            sameSource{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        [perf_diff_same{i}(j,:), ~] = plots.sequentialRemovePredDims(diffSource{i}(j).X_source, ...
            sameSource{i}(j).X_target, diffSource{i}(j).rankRegress.B_, diffSource{i}(j).rankRegress.optDimReducedRankRegress,...
            diffSource{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        [perf_diff_diff{i}(j,:), ~] = plots.sequentialRemovePredDims(diffSource{i}(j).X_source, ...
            diffSource{i}(j).X_target, diffSource{i}(j).rankRegress.B_, diffSource{i}(j).rankRegress.optDimReducedRankRegress,...
            diffSource{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
    end
end

%%
% average and std

for i = 1:nPatterns
    perf_same_same{3+i} = mean(perf_same_same{i},1);
    perf_same_diff{3+i} = mean(perf_same_diff{i},1);
    perf_diff_same{3+i} = mean(perf_diff_same{i},1);
    perf_diff_diff{3+i} = mean(perf_diff_diff{i},1);
    
    perf_same_same{6+i} = std(perf_same_same{i}(1:median_optDim),1);
    perf_same_diff{6+i} = std(perf_same_diff{i}(1:median_optDim),1);
    perf_diff_same{6+i} = std(perf_diff_same{i}(1:median_optDim),1);
    perf_diff_diff{6+i} = std(perf_diff_diff{i}(1:median_optDim),1);
end
%%
% plot
figure (6666)
clf
subplot(3,2,1)
% theta
if numel(perf_same_same{7}) <=1
    plot(perf_same_same{4},'marker','o');
else
    errorbar(perf_same_same{4}, perf_same_same{7},'marker','o');
end
hold on
if numel(perf_same_diff{7}) <=1
    plot(perf_same_diff{4},'marker','o');
else
    errorbar(perf_same_diff{4}, perf_same_diff{7},'marker','o');
end
ylabel("performance")
legend("HPC dims", "PFC dims")
title("predicting target HPC - theta")

subplot(3,2,2)
if numel(perf_diff_same{7}) <=1
    plot(perf_diff_same{4},'marker','o');
else
    errorbar(perf_diff_same{4}, perf_diff_same{7},'marker','o');
end
hold on
if numel(perf_diff_diff{7}) <=1
    plot(perf_diff_diff{4},'marker','o');
else
    errorbar(perf_diff_diff{4}, perf_diff_diff{7},'marker','o');
end
ylabel("performance")
legend("HPC dims", "PFC dims")
title("predicting PFC - theta")

subplot(3,2,3)
% theta
if numel(perf_same_same{8}) <=1
    plot(perf_same_same{5},'marker','o');
else
    errorbar(perf_same_same{5}, perf_same_same{8},'marker','o');
end
hold on
if numel(perf_same_diff{8}) <=1
    plot(perf_same_diff{5},'marker','o');
else
    errorbar(perf_same_diff{5}, perf_same_diff{8},'marker','o');
end
ylabel("performance")
legend("HPC dims", "PFC dims")
title("predicting target HPC - delta")

subplot(3,2,4)
if numel(perf_diff_same{8}) <=1
    plot(perf_diff_same{5},'marker','o');
else
    errorbar(perf_diff_same{5}, perf_diff_same{8},'marker','o');
end
hold on
if numel(perf_diff_diff{8}) <=1
    plot(perf_diff_diff{5},'marker','o');
else
    errorbar(perf_diff_diff{5}, perf_diff_diff{8},'marker','o');
end

ylabel("performance")
legend("HPC dims", "PFC dims")
title("predicting PFC - delta")




subplot(3,2,5)
% ripple
if numel(perf_same_same{9}) <=1
    plot(perf_same_same{6},'marker','o');
else
    errorbar(perf_same_same{6}, perf_same_same{9},'marker','o');
end
hold on
if numel(perf_same_diff{9}) <=1
    plot(perf_same_diff{6},'marker','o');
else
    errorbar(perf_same_diff{6}, perf_same_diff{9},'marker','o');
end
ylabel("performance")
legend("HPC dims", "PFC dims")
title("predicting target HPC - ripple")

subplot(3,2,6)
if numel(perf_diff_same{9}) <=1
    plot(perf_diff_same{6},'marker','o');
else
    errorbar(perf_diff_same{6}, perf_diff_same{9},'marker','o');
end
hold on
if numel(perf_diff_diff{9}) <=1
    plot(perf_diff_diff{6},'marker','o');
else
    errorbar(perf_diff_diff{6}, perf_diff_diff{9},'marker','o');
end

ylabel("performance")
legend("HPC dims", "PFC dims")
title("predicting PFC - ripple")
%%
% pattern-specific subspaces

opt_theta  = cell(1,2);
opt_delta  = cell(1,2);
opt_ripple = cell(1,2);


all_optdim = true; % whether to use a global median for

dimToUse = zeros(Option.numPartition, 2);


for p = 1:Option.numPartition
    for i = 1:2
        if ~all_optdim
            dimToUse(p,i) = max([Patterns(p,i,1).rankRegress.optDimReducedRankRegress,...
                Patterns(p,i,2).rankRegress.optDimReducedRankRegress,...
                Patterns(p,i,3).rankRegress.optDimReducedRankRegress]);
        else
            dimToUse(p,i) = floor(median_optDim);
        end
        
        opt_theta{i}  = [opt_theta{i}, Patterns(p,i,1)];
        opt_delta{i}  = [opt_delta{i}, Patterns(p,i,2)];
        opt_ripple{i} = [opt_ripple{i}, Patterns(p,i,3)];
        
        
    end
end




%%


perf_theta = cell(3,2);
% HPC theta predicting HPC theta; PFC theta predicting PFC theta
% HPC theta predicting HPC delta; PFC theta predicting PFC  delta
% HPC theta predicting HPC ripple; PFC theta predicting PFC ripple
perf_delta = cell(3,2);
perf_ripple = cell(3,2);


for i = 1:2
    for j = 1:nPatterns
        if all_optdim
            perf_theta{j,i} = zeros(numel(opt_theta{i}), floor(median_optDim));
            perf_delta{j,i} = zeros(numel(opt_delta{i}), floor(median_optDim));
            perf_ripple{j,i} = zeros(numel(opt_ripple{i}), floor(median_optDim));
            
        else
            for j = 1:numel(opt_theta)
                % get the largest optDim for theta?
            end
            perf_theta = zeros(numel(opt_theta{i}), max(opt_theta{i}));
            perf_delta = zeros(numel(opt_delta{i}), floor(median_to_use(i)));
            perf_ripple = zeros(numel(opt_ripple{i}), floor(median_to_use(i)));
        end
    end
end

%%

for i = 1:2
    
    
    for j = 1:num_same
        
        % theta-theta, delta-delta, ripple-ripple
        [perf_theta{1,i}(j,:), ~] = plots.sequentialRemovePredDims(opt_theta{i}(j).X_source, ...
            opt_theta{i}(j).X_target, opt_theta{i}(j).rankRegress.B_, median_optDim,...
            opt_theta{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        
        [perf_delta{2,i}(j,:), ~] = plots.sequentialRemovePredDims(opt_delta{i}(j).X_source, ...
            opt_delta{i}(j).X_target, opt_delta{i}(j).rankRegress.B_, median_optDim,...
            opt_delta{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        
        [perf_ripple{3,i}(j,:), ~] = plots.sequentialRemovePredDims(opt_ripple{i}(j).X_source, ...
            opt_ripple{i}(j).X_target, opt_ripple{i}(j).rankRegress.B_, median_optDim,...
            opt_ripple{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        
        % theta-delta, delta-theta, ripple-theta
        [perf_theta{2,i}(j,:), ~] = plots.sequentialRemovePredDims(opt_delta{i}(j).X_source, ...
            opt_delta{i}(j).X_target, opt_theta{i}(j).rankRegress.B_, median_optDim,...
            opt_theta{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        
        [perf_delta{1,i}(j,:), ~] = plots.sequentialRemovePredDims(opt_theta{i}(j).X_source, ...
            opt_theta{i}(j).X_target, opt_delta{i}(j).rankRegress.B_, median_optDim,...
            opt_delta{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        
        [perf_ripple{1,i}(j,:), ~] = plots.sequentialRemovePredDims(opt_theta{i}(j).X_source, ...
            opt_theta{i}(j).X_target, opt_ripple{i}(j).rankRegress.B_, median_optDim,...
            opt_ripple{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        
        % theta-ripple, delta-ripple, ripple-delta
        [perf_theta{3,i}(j,:), ~] = plots.sequentialRemovePredDims(opt_ripple{i}(j).X_source, ...
            opt_ripple{i}(j).X_target, opt_theta{i}(j).rankRegress.B_, median_optDim,...
            opt_theta{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        
        [perf_delta{3,i}(j,:), ~] = plots.sequentialRemovePredDims(opt_ripple{i}(j).X_source, ...
            opt_ripple{i}(j).X_target, opt_delta{i}(j).rankRegress.B_, median_optDim,...
            opt_delta{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        
        [perf_ripple{2,i}(j,:), ~] = plots.sequentialRemovePredDims(opt_delta{i}(j).X_source, ...
            opt_delta{i}(j).X_target, opt_ripple{i}(j).rankRegress.B_, median_optDim,...
            opt_ripple{i}(j).rankRegress.cvLoss, numDimsUsedForPrediction, "normalized", false);
        
    end
end

%%
% average and std
for p = 1:3
    for i = 1:2
        perf_theta{p,i}(Option.numPartition+1,:) = mean(perf_theta{p,i}(1:median_optDim,:),1);
        perf_delta{p,i}(Option.numPartition+1,:) = mean(perf_delta{p,i}(1:median_optDim,:),1);
        perf_ripple{p,i} (Option.numPartition+1,:)= mean(perf_ripple{p,i}(1:median_optDim,:),1);
        
        
        perf_theta{p,i} (Option.numPartition+2,:)= std(perf_theta{p,i}(1:median_optDim,:),1);
        perf_delta{p,i} (Option.numPartition+2,:)= std(perf_delta{p,i}(1:median_optDim,:),1);
        perf_ripple{p,i} (Option.numPartition+2,:)= std(perf_ripple{p,i}(1:median_optDim,:),1);
        
    end
end

%%
figure (6667)
clf


subplot(3,2,1)
% theta - HPC: predicting its theta, delta, and ripple
errorbar(perf_theta{1,1}(Option.numPartition+1,:), perf_theta{1,1}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_theta{2,1}(Option.numPartition+1,:), perf_theta{2,1}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_theta{3,1}(Option.numPartition+1,:), perf_theta{3,1}(Option.numPartition+2,:),'marker','o');

ylabel("performance")
legend("theta dims", "delta dims", "ripple dims")
title("predicting target HPC using theta")

subplot(3,2,2)
% theta - HPC: predicting its theta, delta, and ripple
errorbar(perf_theta{1,2}(Option.numPartition+1,:), perf_theta{1,2}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_theta{2,2}(Option.numPartition+1,:), perf_theta{2,2}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_theta{3,2}(Option.numPartition+1,:), perf_theta{3,2}(Option.numPartition+2,:),'marker','o');

ylabel("performance")
legend("theta dims", "delta dims", "ripple dims")
title("predicting target PFC using theta")

subplot(3,2,3)
errorbar(perf_delta{1,1}(Option.numPartition+1,:), perf_delta{1,1}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_delta{2,1}(Option.numPartition+1,:), perf_delta{2,1}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_delta{3,1}(Option.numPartition+1,:), perf_delta{3,1}(Option.numPartition+2,:),'marker','o');

ylabel("performance")
legend("theta dims", "delta dims", "ripple dims")
title("predicting target HPC using delta")

subplot(3,2,4)
errorbar(perf_delta{1,2}(Option.numPartition+1,:), perf_delta{1,2}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_delta{2,2}(Option.numPartition+1,:), perf_delta{2,2}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_delta{3,2}(Option.numPartition+1,:), perf_delta{3,2}(Option.numPartition+2,:),'marker','o');

ylabel("performance")
legend("theta dims", "delta dims", "ripple dims")
title("predicting target PFC using delta")


subplot(3,2,5)
errorbar(perf_ripple{1,1}(Option.numPartition+1,:), perf_ripple{1,1}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_ripple{2,1}(Option.numPartition+1,:), perf_ripple{2,1}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_ripple{3,1}(Option.numPartition+1,:), perf_ripple{3,1}(Option.numPartition+2,:),'marker','o');

ylabel("performance")
legend("theta dims", "delta dims", "ripple dims")
title("predicting target HPC using ripple")




subplot(3,2,6)
errorbar(perf_ripple{1,2}(Option.numPartition+1,:), perf_ripple{1,2}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_ripple{2,2}(Option.numPartition+1,:), perf_ripple{2,2}(Option.numPartition+2,:),'marker','o');
hold on
errorbar(perf_ripple{3,2}(Option.numPartition+1,:), perf_ripple{3,2}(Option.numPartition+2,:),'marker','o');

ylabel("performance")
legend("theta dims", "delta dims", "ripple dims")
title("predicting target PFC using ripple")

set(findobj(gcf,'type','axes'), 'yscale', 'log')
%%