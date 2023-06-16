%fig('Mean strength stem plot');clf; tiledlayout(2,2);
titles = ["reward", "error", "inBoundChoice","outBoundChoice"];
fig('instaneous correlation of learning and component strength - inbound'); clf;tiledlayout(2,2);
for i = 1:4
    nexttile
    if i == 4
        continue
    end
    for j = 1:30
        hold on
        stem(j:j,critical_behaviors(i).inBound_stable_corr(j),"filled", "color", "black");
        hold on
        stem(j:j,critical_behaviors(i).inBound_learning_corr(j),"color", "black");
        title(titles(i));
        
    end
    hold on
    plots.equalPatches
    hold on
    legend('stable','learning')
    title(titles(i));
    xticks([3,8,13,18,23,28])
    xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
    ylabel("correlation with learning")
end
%%
fig('instaneous correlation of learning and component strength - outbound'); clf
for i = 1:4
    nexttile
    if i == 3
        continue
    end
    for j = 1:30
        hold on
        stem(j:j,critical_behaviors(i).outBound_stable_corr(j),"filled", "color", "black");
        hold on
        stem(j:j,critical_behaviors(i).outBound_learning_corr(j),"color", "black");
        title(titles(i));
        
    end
    hold on
    plots.equalPatches
    hold on
    legend('stable','learning')
    title(titles(i));
    xticks([3,8,13,18,23,28])
    xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
    ylabel("correlation with learning")
end