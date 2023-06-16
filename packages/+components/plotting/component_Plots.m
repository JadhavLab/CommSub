
n = repmat({1:Option.dimCompAnalysis*numel(patternNames)}, 1 ,4);
norm = @(x) (x-nanmin(x))./(nanmax(x)-nanmin(x));
%% first, check if events occur more often during certain ocsillation patterns
%norm no longer works??
fig('Mean strength stem plot'); clf; tiledlayout(2,2);
for i = 1:4
    Z{i} = critical_behaviors(i).count_times;
end

titles = ["error", "inBoundChoice","outBoundChoice","reward"];
for i = 1:4
    nexttile;
    stem(1:7, critical_behaviors(i).count_times);
    xticks(1:7);
    xticklabels(["none", "theta", "delta","ripple", "low-theta","low-delta","low-ripple"])
    title(titles(i));
    
end
%% compare the mean component strength for different rhythms (of a selected dimension)

fig('Mean strength stem plot');clf; tiledlayout(2,2);
titles = ["reward", "error", "inBoundChoice","outBoundChoice"];

for i = 1:4
    nexttile
    for j = 1:5
        hold on
        stem(1:6, critical_behaviors(i).mean_component_strength(j,:));
        xticks(1:6);
        xticklabels(["theta", "delta","ripple", "low-theta","low-delta","low-ripple"])
        title(titles(i));
        
    end
    legend("dim1", "dim2", "dim3", "dim4", "dim5")
end

%%
plot_component_strength

%%
f=fig('Strength of all subspaces with inbound tperf, by behavior ' + animal); clf; tiledlayout(2,2);
titles = ["reward", "error", "inBoundChoice","outBoundChoice"];
ax = gobjects(4);
for i = 1:4
    % restructure components using times
    %------------------------------------
    % get times with jump
%     possible_jumps = find( [0; diff(critical_times{i})] > 60*5 );
%     possible_jumps = max([possible_jumps, possible_jumps-1], 0);
    critical_components = Patterns(10,2,1).compAnalysisByBeh(i).critical_components;
    %critical_components(possible_jumps(:),:) = nan;
    critical_times = Patterns(10,2,1).compAnalysisByBeh(i).critical_times;
    new_time = min(critical_times):(min(diff(critical_times))):max(critical_times);
    for j = 1:size(critical_components,2)
        cc(:,j) = interp1(critical_times, critical_components(:,j), new_time);
    end
    
    % plot
    % ----
    nexttile;
     X = (cc' - nanmedian(cc',2))./nanstd(cc',0,2);
     X = components.organizeComponents(X, 5, 3, 'strength');
     imagesc(new_time, n{i}, X)
    crameri bam
    clim = get(gca,'clim');
    hold on
    if ~isempty(critical_behaviors(i).tperfInbound)
        plot(critical_times, 30*norm(Patterns(10,2,1).compAnalysisByBeh(i).unifiedTperfInbound),'LineWidth',4)
       
    end
    hold on
    if ~isempty(critical_behaviors(i).tperfOutbound)
        plot(critical_times, 30*norm(Patterns(10,2,1).compAnalysisByBeh(i).unifiedTperfOutbound),'LineWidth',4)
        
    end
    legend("inbound","outbound")
    X_time_mean_abs = mean(abs(X),1);
%     box(X_time_mean_abs, 
%     scatter(new_time, 90*norm(X_time_mean_abs), 8, 'MarkerFaceColor','b','MarkerEdgeColor','b',...
%            'MarkerFaceAlpha',.01,'MarkerEdgeAlpha',.01);
   % legend("inbound perf", "outbound perf")
        nanloc = isnan(cc(:,1));
%     a=area(new_time, 31*nanloc, 'FaceColor', 'black', 'FaceAlpha', 0.8);
    
    set(gca,'ydir','normal')
    ax(i) = gca;
  
    hold on
    title(titles(i));
      xlabel("time")
    ylabel("component")
    yticks([3,8,13,18,23,28])
    yticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
%     yticks([1,6,11,16,21,26])
%     yticklabels(["theta", "delta","ripple", "low-theta","low-delta","low-ripple"])
end
% linkaxes(ax, 'xy')

%% correlation of all the components, by critical behaviors
plotCpCorr
%%
plotCpStrengthLearning
%%
plotInstantLearningCorr