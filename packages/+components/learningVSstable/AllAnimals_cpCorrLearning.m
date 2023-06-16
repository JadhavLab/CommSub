% after calculating the subspace component over all animals

% CURRENTLY UNAVAILABLE UNTIL CORR IS FIGURED OUT

%% looks at the sum of the learning vs. stable component strengths for all animals

for i = 1:numel(animal_list)
    for j = 1:4
        
        mean_inbound_learning_component{i}{j} = nanmean(animal_critical_behaviors{i}(j).inBound_learning_component,1);
        mean_outbound_learning_component{i}{j}= nanmean(animal_critical_behaviors{i}(j).outBound_learning_component,1);
        
        mean_inbound_stable_component{i}{j} = nanmean(animal_critical_behaviors{i}(j).inBound_stable_component,1);
        mean_outbound_stable_component{i}{j} = nanmean(animal_critical_behaviors{i}(j).outBound_stable_component,1);
        
        mean_inbound_learning(i,j) = mean(mean_inbound_learning_component{i}{j});
        mean_inbound_stable(i,j) = mean(mean_inbound_stable_component{i}{j} );
        
        mean_outbound_learning(i,j) = mean(mean_outbound_learning_component{i}{j});
        mean_outbound_stable(i,j) = mean(mean_outbound_stable_component{i}{j});
        
    end
end

%% plot
fig('mean component strength, all patterns, during different behaviors'); clf; tiledlayout(2,2);
inbound_learning = nanmean(abs(mean_inbound_learning),1);
inbound_stable = nanmean(abs(mean_inbound_stable),1);
outbound_learning = nanmean(abs(mean_outbound_learning),1);
outbound_stable = nanmean(abs(mean_outbound_stable),1);


std_inbound_learning = nanstd(abs(mean_inbound_learning),1);
std_inbound_stable = nanstd(abs(mean_inbound_stable),1);
std_outbound_learning = nanstd(abs(mean_outbound_learning),1);
std_outbound_stable = nanstd(abs(mean_outbound_stable),1);
titles = ["reward", "error", "inBoundChoice","outBoundChoice"];

for i = 1:4
    nexttile;
    bar([inbound_learning(i),inbound_stable(i),outbound_learning(i),outbound_stable(i)]);
    hold on
    errorbar([inbound_learning(i),inbound_stable(i),outbound_learning(i),outbound_stable(i)],...
        [std_inbound_learning(i),std_inbound_stable(i),std_outbound_learning(i),std_outbound_stable(i)]);
    
    xticks([1,2,3,4])
    xticklabels(["inbound learning", "inbound stable","outbound learning", "outbound stable"])
    title(titles(i));
    ylabel("average component strength")
    
end
%%
fig('mean component strength, by patterns, during different behaviors of learning vs stable, inbound');
clf; tiledlayout(2,2);

for i = 1:4
    
    nexttile
    for k = 1:numel(animal_list)
        if i == 4
            continue
        end
        for j = 1:30
            hold on
            stem(j:j, mean_inbound_stable_component{k}{i}(j) ,"filled", "color", "black");
            hold on
            stem(j:j,mean_inbound_learning_component{k}{i}(j),"color", "black");
            
            
        end
     
        
    end
       hold on
        plots.equalPatches
        hold on
        legend('stable','learning')
        title(titles(i));
        xticks([3,8,13,18,23,28])
        xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
        ylabel("mean component strength")
    title(titles(i));
end
%% just plot the difference to see variability...
fig('mean component strength, by patterns, during different behaviors of learning vs stable, inbound');
clf; tiledlayout(2,2);

for i = 1:4
    
    nexttile
    for k = 1:numel(animal_list)
        if i == 4
            continue
        end
        for j = 1:30

            stem(j:j, abs(mean_inbound_stable_component{k}{i}(j))- ...
                abs(mean_inbound_learning_component{k}{i}(j)),"filled", "color", "black");
            hold on

        end
    end
       hold on
        plots.equalPatches
        hold on

        title(titles(i));
        xticks([3,8,13,18,23,28])
        xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
        ylabel("mean component strength")
    title(titles(i));
end


%%
fig('mean component strength, by patterns, during different behaviors of learning vs stable, outbound');
clf; tiledlayout(2,2);

for i = 1:4
    
    nexttile
    for k = 1:numel(animal_list)
        if i == 3
            continue
        end
        for j = 1:30
            hold on
            stem(j:j, mean_outbound_stable_component{k}{i}(j) ,"filled", "color", "black");
            hold on
            stem(j:j,mean_outbound_learning_component{k}{i}(j),"color", "black");
            
            
        end
     
        
    end
       hold on
        plots.equalPatches
        hold on
        legend('stable','learning')
        title(titles(i));
        xticks([3,8,13,18,23,28])
        xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
        ylabel("mean component strength")
    title(titles(i));
end


%
%
%             mean_outbound_learning_component{i}{j}
%             mean_outbound_stable_component{i}{j}
%%
fig('inst corr, by patterns, during different behaviors of learning vs stable, inbound');
clf; tiledlayout(2,2);

for i = 1:4
    
    nexttile
    for k = 1:numel(animal_list)
        if i == 4
            continue
        end
        for j = 1:30
            hold on
%             temp_max = max(max(animal_critical_behaviors{k}(i).outBound_stable_corr),...
%                 max(animal_critical_behaviors{k}(i).outBound_learning_corr));
            stem(j:j,abs(animal_critical_behaviors{k}(i).inBound_stable_corr(j))-...
                abs(animal_critical_behaviors{k}(i).inBound_learning_corr(j)),"filled","color", "black");
        end
    end
       hold on
        plots.equalPatches
        hold on
        title(titles(i));
        xticks([3,8,13,18,23,28])
        xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
        ylabel("inst corr (stable - learning)")
    title(titles(i));
end
%%
fig('inst corr, by patterns, during different behaviors of learning vs stable, outbound');
clf; tiledlayout(2,2);

for i = 1:4
    
    nexttile
    for k = 1:numel(animal_list)
        if i == 3
            continue
        end
        for j = 1:30
            hold on
%             temp_max = max(max(animal_critical_behaviors{k}(i).outBound_stable_corr),...
%                 max(animal_critical_behaviors{k}(i).outBound_learning_corr));
            stem(j:j,abs(animal_critical_behaviors{k}(i).outBound_stable_corr(j))-...
                abs(animal_critical_behaviors{k}(i).outBound_learning_corr(j)),"filled","color", "black");

        end
     
        
    end
       hold on
        plots.equalPatches
        hold on
        title(titles(i));
        xticks([3,8,13,18,23,28])
        xticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
        ylabel("inst corr (stable - learning)")
    title(titles(i));
end


%%


