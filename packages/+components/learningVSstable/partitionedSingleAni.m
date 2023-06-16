titles = ["reward", "error", "choice"];
for p = 1:Option.numPartition
    hpcCompBeh = Patterns(p,1,1).compAnalysisByBeh;
    pfcCompBeh = Patterns(p,2,1).compAnalysisByBeh;
    
    
    for i = 1:4
        if i == 4
            continue;
        end
        hpcInDiff{i}(p,:) = abs(hpcCompBeh(i).inBound_stable_corr) - ...
            abs(hpcCompBeh(i).inBound_learning_corr);
        pfcInDiff{i}(p,:) = abs(pfcCompBeh(i).inBound_stable_corr) - ...
            abs(pfcCompBeh(i).inBound_learning_corr);
    end
    
end

for p = 1:Option.numPartition
    hpcCompBeh = Patterns(p,1,1).compAnalysisByBeh;
    pfcCompBeh = Patterns(p,2,1).compAnalysisByBeh;
    
    
    for i = 1:4
        if i == 3
            continue;
        end
        hpcOutDiff{i}(p,:) = abs(hpcCompBeh(i).outBound_stable_corr) - ...
            abs(hpcCompBeh(i).outBound_learning_corr);
        pfcOutDiff{i}(p,:) = abs(pfcCompBeh(i).outBound_stable_corr) - ...
            abs(pfcCompBeh(i).outBound_learning_corr);
    end
    
end
%%
for j = 1:4
    if j == 4
        mean_hpcOutDiff(j,:) = mean(hpcOutDiff{j},1);
        mean_pfcOutDiff(j,:) = mean(pfcOutDiff{j},1);
    elseif j ==3
        mean_hpcInDiff(j,:) = mean(hpcInDiff{j},1);
        mean_pfcInDiff(j,:) = mean(pfcInDiff{j},1);
    else
        mean_hpcInDiff(j,:) = mean(hpcInDiff{j},1);
        mean_pfcInDiff(j,:) = mean(pfcInDiff{j},1);
        mean_hpcOutDiff(j,:) = mean(hpcOutDiff{j},1);
        mean_pfcOutDiff(j,:) = mean(pfcOutDiff{j},1);
    end
end

%%
figure
fig('inst corr, by patterns, during different behaviors of learning vs stable, inbound');
clf; tiledlayout(2,2);

for i = 1:4
    
    nexttile
    
    if i == 4
        continue
    end
    for j = 1:30
        hold on
        %             temp_max = max(max(animal_critical_behaviors{k}(i).outBound_stable_corr),...
        %                 max(animal_critical_behaviors{k}(i).outBound_learning_corr));
        stem(j:j,mean_hpcInDiff(i,j),"filled","color", "black");
        hold on
        stem(j:j,mean_pfcInDiff(i,j),"color", "black");
        
    end
    
    hold on
    plots.equalPatches
    hold on
    legend("hpc-hpc","hpc-pfc")
    title(titles(i));
    xticks([3,8,13,18,23,28])
    xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
    ylabel("inst corr (stable - learning)")
    
end


%%
fig('inst corr, by patterns, during different behaviors of learning vs stable, outbound');
clf; tiledlayout(2,2);

for i = 1:4
    
    nexttile
    
    if i == 3
        continue
    end
    for j = 1:30
        hold on
        %             temp_max = max(max(animal_critical_behaviors{k}(i).outBound_stable_corr),...
        %                 max(animal_critical_behaviors{k}(i).outBound_learning_corr));
        stem(j:j,mean_hpcOutDiff(i,j),"filled","color", "black");
        hold on
        stem(j:j,mean_pfcOutDiff(i,j),"color", "black");
        
    end
    
    hold on
    plots.equalPatches
    hold on
    legend("hpc-hpc","hpc-pfc")
    title(titles(i));
    xticks([3,8,13,18,23,28])
    xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
    ylabel("inst corr (stable - learning)")
    
end