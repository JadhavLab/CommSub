% give Pattern struct, compute the #dims spanned by all the partitions of
% different rhythms

% Patterns.name

%% sum all, not particularly meaningful.
patternDimSum = zeros(1,nPatterns);
for a = 1:nAnimal
    for m = 1:nMethods
        for p = 1:nPartition
            for d = 1:2
                for i = 1:nPatterns
                    curr = Patterns(a,m,p,d,i);
                    if curr.name == "theta"
                        patternDimSum(1) = patternDimSum(1) ...
                            + curr.rankRegress.optDimReducedRankRegress;
                    elseif curr.name == "delta"
                        patternDimSum(2) = patternDimSum(2) ...
                            + curr.rankRegress.optDimReducedRankRegress;
                    elseif curr.name == "ripple"
                        patternDimSum(3) = patternDimSum(3) ...
                            + curr.rankRegress.optDimReducedRankRegress;
                    end
                end
            end
        end
    end
end

%% percent by patterns - significant among rhythms
patternPercSum = zeros(nAnimal,2,nPatterns);
all_theta = [];
all_delta = [];
all_ripple = [];

for a = 1:nAnimal
    for m = 1:1
        for p = 1:nPartition
            for d = 1:2
                for i = 1:nPatterns
                    curr = Patterns(a,m,p,d,i);
                    nTarget_temp = numel(curr.index_target);
                    nSource_temp = numel(curr.index_source);
                    nDimMax = min(nTarget_temp, nSource_temp);
                    currDim = curr.rankRegress.optDimReducedRankRegress;
                    currPerc = currDim/nDimMax;
                    if curr.name == "theta"
                        all_theta = [all_theta, currPerc];
                        patternPercSum(a,d,1) = patternPercSum(a,d,1) ...
                            + currPerc;
                    elseif curr.name == "delta"
                        all_delta = [all_delta, currPerc];
                        patternPercSum(a,d,2) = patternPercSum(a,d,2) ...
                            + currPerc;
                    elseif curr.name == "ripple"
                        all_ripple = [all_ripple, currPerc];
                        patternPercSum(a,d,3) = patternPercSum(a,d,3) ...
                            + currPerc;
                    end
                end
            end
        end
    end
end

%% plot

fig("comparing rhythm dimensions")

% bar plots of all rhythm mean across animals and std
x = 1:3;
mean_patternPerc = squeeze(mean(patternPercSum,1));
std_patternPerc = squeeze(std(patternPercSum,1));
bar(x,[mean_patternPerc(1,:); mean_patternPerc(2,:)])
alpha(0.33)
hold on
er = errorbar(x,mean_patternPerc(1,:), std_patternPerc(1,:));
er2 = errorbar(x,mean_patternPerc(2,:), std_patternPerc(2,:));

xticklabels(patternnames)
ylabel("sum of perc. dims spanned")

legend("hpc-hpc","hpc-pfc")


[h_compare, p_compare] = ttest2(mean_patternPerc(1,:), mean_patternPerc(2,:))

% dim perc across patterns
temp = permute(patternPercSum,[2,1,3]);
temp_hpc = temp(1,:,:);
temp_pfc = temp(2,:,:);
hpc_percdim = squeeze(temp_hpc);
pfc_percdim = squeeze(temp_pfc);

[p_hpc] = anova1(hpc_percdim)
[p_pfc] = anova1(pfc_percdim)
%% percent by actual/control
timePercSum = zeros(nAnimal,2,nPatterns); % the two here is high time/control

all_control = cell(1,3);
all_acti = cell(1,3);

for i = 1:3
    all_control{i} = [];
    all_acti{i} = [];
end

for a = 1:nAnimal
    for m = 2:2
        for p = 1:nPartition
            for d = 2:2
                for i = 1:nPatterns*2
                    curr = Patterns(a,m,p,d,i);
                    nTarget_temp = numel(curr.index_target);
                    nSource_temp = numel(curr.index_source);
                    nDimMax = min(nTarget_temp, nSource_temp);
                    currDim = curr.rankRegress.optDimReducedRankRegress;
                    currPerc = currDim/nDimMax;
                    if i/3 >1 % control
                        if curr.name == "theta-control"
                            all_control{1} = [all_control{1}, currPerc];
                            timePercSum(a,2,1) = timePercSum(a,2,1) ...
                                + currPerc;
                        elseif curr.name == "ripple-control"
                            all_control{2} = [all_control{2}, currPerc];
                            timePercSum(a,2,2) = timePercSum(a,2,2) ...
                                + currPerc;
                        elseif curr.name == "delta-control"
                            all_control{3} = [all_control{3}, currPerc];
                            timePercSum(a,2,3) = timePercSum(a,2,3) ...
                                + currPerc;
                        end
                    else
                        if curr.name == "theta"
                            all_acti{1}  = [all_acti{1} , currPerc];
                            timePercSum(a,1,1) = timePercSum(a,1,1) ...
                                + currPerc;
                        elseif curr.name == "delta"
                            all_acti{2}  = [all_acti{2} , currPerc];
                            timePercSum(a,1,2) = timePercSum(a,1,2) ...
                                + currPerc;
                        elseif curr.name == "ripple"
                            all_acti{3}  = [all_acti{3} , currPerc];
                            timePercSum(a,1,3) = timePercSum(a,1,3) ...
                                + currPerc;
                        end
                    end
                end
            end
        end
    end
end

%% plot

fig("comparing activity dimensions"), clf

% bar plots of all rhythm mean across animals and std
x = 1:3;
mean_timePerc = squeeze(mean(timePercSum,1));
std_timePerc = squeeze(std(timePercSum,1));
bar(x,[mean_timePerc(1,:); mean_timePerc(2,:)])
alpha(0.33)
hold on
er = errorbar(x,mean_timePerc(1,:), std_timePerc(1,:));
er2 = errorbar(x,mean_timePerc(2,:), std_timePerc(2,:));

xticklabels(patternnames)
ylabel("sum of perc. dims spanned")

legend("activity","control")


control = [];
acti = [];
for i = 1:3
    [h_pattern(i), p_pattern(i)] = ttest2(all_control{i}, all_acti{i});
    control = [control, all_control{i}];
    acti = [acti, all_acti{i}];
    [h_compare(i), p_compare(i)] = ttest2(control, acti)
end

