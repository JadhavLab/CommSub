%%
%mean component strength during
for i = 1:4
    mean_inbound_learning_component{i} = nanmean(critical_behaviors(i).inBound_learning_component,1);
    mean_outbound_learning_component{i} = nanmean(critical_behaviors(i).outBound_learning_component,1);
    
    mean_inbound_stable_component{i} = nanmean(critical_behaviors(i).inBound_stable_component,1);
    mean_outbound_stable_component{i} = nanmean(critical_behaviors(i).outBound_stable_component,1);
    
end

fig('mean component strength during learning and stable phase - inbound'); 
clf; tiledlayout(2,2);
for i = 1:4
    
    nexttile
    if i == 4
        continue
    end
    for j = 1:30
        hold on
        stem(j:j,mean_inbound_stable_component{i}(j),"filled", "color", "black");
        hold on
        stem(j:j,mean_inbound_learning_component{i}(j),"color", "black");
        title(titles(i));
        
    end
    hold on
    plots.equalPatches
    hold on
    legend('learning','stable')
    title(titles(i));
    xticks([3,8,13,18,23,28])
    xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
    ylabel("mean component strength")
    
end


fig('mean component strength during learning and stable phase - outbound'); 
clf;tiledlayout(2,2);
for i = 1:4
    nexttile
    if i == 3
        continue
    end
    
    
    for j = 1:30
        hold on
        stem(j:j,mean_outbound_stable_component{i}(j),"filled", "color", "black");
        hold on
        stem(j:j,mean_outbound_learning_component{i}(j),"color", "black");
        title(titles(i));
        
    end
    hold on
    plots.equalPatches
    hold on
    legend('learning','stable')
    title(titles(i));
    xticks([3,8,13,18,23,28])
    xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
    ylabel("mean component strength")
    
end
%%
for i = 1:4
    mean_inbound_learning(i) = sum(mean_inbound_learning_component{i});
    mean_inbound_stable(i) = sum(mean_outbound_stable_component{i});
    
    mean_outbound_learning(i) = sum(mean_outbound_learning_component{i});
    mean_outbounr_stable(i) = sum(mean_outbound_stable_component{i});
end