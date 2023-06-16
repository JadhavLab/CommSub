%% unified time with component strengths across time


titles = ["reward", "error", "inBoundChoice","outBoundChoice"];
%
fig('Temporal presence of factors, all times');
clf;
tiledlayout(2,2);
for i = 1:4
    
    nexttile; imagesc(Patterns(1,2,1).compAnalysisByBeh(i).critical_times, n{i}, Patterns(1,1,1).compAnalysisByBeh(i).critical_components')
    
    for j = 1:5
        lineObject = line([0 max(critical_times{i})] ,[j*5 j*5]);
        lineObject.LineStyle = ':'; % Make line dotted
        lineObject.LineWidth = 2;  % Thicken the line
        lineObject.Color = 'white'; %
        hold on
    end
    title(titles(i));
    xlabel("time")
    ylabel("component")
end

%% Strength of all subspaces during different critical time points
% RY - median scaled, and then clim centeterd
fig('Temporal presence of factors, median centered, each critical times');clf;
tiledlayout(2,2);
for i = 1:4
    Y_original{i} = Patterns(1,1,1).compAnalysisByBeh(i).critical_components;
end
titles = ["reward", "error", "inBoundChoice","outBoundChoice"];
sgtitle('X - median(X) and centered color scale')

for i = 1:4
    nexttile;
    imagesc(critical_times{i}, n{i}, Y_original{i}' - nanmedian(Y_original{i}',2))
    %cmocean haline
    crameri bam
    clim = get(gca,'clim');
    
    for j = 1:5
        lineObject = line ([0 max(critical_times{i})] ,[j*5+0.5 j*5+0.5]);
        lineObject.LineStyle = ':'; % Make line dotted
        lineObject.LineWidth = 2;  % Thicken the line
        lineObject.Color = 'black'; %
        hold on
    end
    title(titles(i));
    xlabel("time")
    ylabel("component")
end
set(findobj(gcf,'type','axes'), 'clim', [-max(clim(:)), max(clim(:))])
sgtitle('X - median(X) and centered color scale')


%% Strength of all subspaces during different critical time points
% RY - median scaled, and then clim centeterd
fig('Temporal presence of factors, median centered, std scaled, each critical time'); clf; tiledlayout(2,2);

titles = ["reward", "error", "inBoundChoice","outBoundChoice"];

criticalTimes = {animal_behavior.rewardTimes, animal_behavior.errorTimes,...
    animal_behavior.inBoundChoiceTimes, animal_behavior.outBoundChoiceTimes};

for i = 1:4
    Y_original{i} = components.organizeComponents(Y_original{i}', 5, 3, "strength");
end
for i = 1:4
    nexttile;
    imagesc(critical_times{i}, n{i}, (Y_original{i} - nanmedian(Y_original{i},2))./std(Y_original{i},0,2))
    %cmocean haline
    crameri bam
    clim = get(gca,'clim');
    
    for j = 1:5
        lineObject = line ([0 max(critical_times{i})] ,[j*5+0.5 j*5+0.5]);
        lineObject.LineStyle = ':'; % Make line dotted
        lineObject.LineWidth = 2;  % Thicken the line
        lineObject.Color = 'black'; %
        hold on
    end
    title(titles(i));
    xlabel("time")
    ylabel("component")
    yticks([1,6,11,16,21,26])
    yticklabels(["theta", "low-theta","delta","low-delta","ripple","low-ripple"])
end
set(findobj(gcf,'type','axes'), 'clim', [-max(clim(:)), max(clim(:))])
sgtitle('X - median(X) and centered color scale, divided by std of each pat component')

%% Strength of all subspaces with lindist plotted on top
fig('Strength of all subspaces with lindist' + animal); clf;
% set(gcf,'position', get(0, 'screensize'));
imagesc(Patterns(1,2,1).allSubSpaceComponents)
crameri bam
clim = get(gca,'clim');
hold on

lindist_to_plot = ((animal_behavior.lindist-median(animal_behavior.lindist))...
    /(max(animal_behavior.lindist)-median(animal_behavior.lindist)))' * 30;
plot(lindist_to_plot, '-k')


%% Strength of all subspaces with tperf (in/outBound) plotted on top

fig('Strength of all subspaces with inbound tperf'); clf;
% Rset(gcf,'position', get(0, 'screensize'));
imagesc(Patterns(1,2,1).allSubSpaceComponents)
crameri bam
clim = get(gca,'clim');
hold on

norm = @(x) (x-nanmin(x))./(nanmax(x)-nanmin(x));
inBoundPerf_to_plot = 30*(norm(animal_behavior.tperfInbound));
plot(inBoundPerf_to_plot, '-k')
set(gca,'ydir','normal')
hold on
outBoundPerf_to_plot = 30*(norm(animal_behavior.tperfOutbound));
plot(outBoundPerf_to_plot, '-b')
xlabel("time")
    yticks([1,6,11,16,21,26])
    yticklabels(["theta", "delta","ripple","low-theta","low-delta","low-ripple"])
set(gca,'ydir','normal')
legend('inbound', 'outbound')
