function [critical_behaviors, unifiedTime, critical_components, critical_times] = makeBehaviorStructs...
    (dimension, cellOfWindows, animal_behavior,  behavior_running_subspaces, unify_from_running)

% this function makes structs for the four types of critical behaviors
% interested in

% input:
% dimension: the dimensions of a single pattern analyzed
% cellOfWindows: for analyzing the actual event happening
% animal_behavior: behavior lookup table for the given animal
% unify_from_running: whether the unified time is all that when the animal
% is running, or just unifying all the critical behavioral times
% quantile: the quantile for determining whether the animal is learning or
% has achieved stable performance

% output:
% critical_behaviors: an array of structs, each containing variables of interests regrading
% a particular behavior of interest
% unifiedTime: the unified time for all behavior of interest


rewarded = struct();
errored = struct();
inBoundChoice = struct();
outBoundChoice = struct();

critical_behaviors = [rewarded, errored, inBoundChoice, outBoundChoice];
subspaceComponents = behavior_running_subspaces';
%% where events occur

% sample all the time points and find the type of event they belong,
% eventTypes = components.addEventType(cellOfWindows, animal_behavior);
% animal_behavior.eventTypes = eventTypes';
% count_inBoundChoiceTimes = zeros(1,7);
% count_outBoundChoiceTimes = zeros(1,7);
% count_errorTimes = zeros(1,7);
% count_rewardTimes = zeros(1,7);
% 
% 
% % count the number of times special events occur during each pattern
% % activity
% for i = 1:7
%     count_inBoundChoiceTimes(i) = sum(animal_behavior.inBoundChoiceTimes & ...
%         animal_behavior.eventTypes ==i);
%     count_outBoundChoiceTimes(i) = sum(animal_behavior.outBoundChoiceTimes & ...
%         animal_behavior.eventTypes ==i);
%     count_errorTimes(i) = sum(animal_behavior.errorTimes & ...
%         animal_behavior.eventTypes ==i);
%     count_rewardTimes(i) = sum(animal_behavior.rewardTimes & ...
%         animal_behavior.eventTypes ==i);
% end
% 
% K = {count_rewardTimes,count_errorTimes,count_inBoundChoiceTimes,count_outBoundChoiceTimes};
% for i = 1:4
%     critical_behaviors(i).count_times = K{i};
% end

%% component strength for different rhythms
rewardedComponent           = subspaceComponents(animal_behavior.rewardTimes,:);
erroredComponent            = subspaceComponents(animal_behavior.errorTimes,:);
inBoundChoiceComponent      = subspaceComponents(animal_behavior.inBoundChoiceTimes,:);
outBoundChoiceComponent     = subspaceComponents(animal_behavior.outBoundChoiceTimes,:);

rewardedComponentTime       = animal_behavior.time(animal_behavior.rewardTimes);
erroredComponentTime        = animal_behavior.time(animal_behavior.errorTimes);
inBoundChoiceComponentTime  = animal_behavior.time(animal_behavior.inBoundChoiceTimes);
outBoundChoiceComponentTime = animal_behavior.time(animal_behavior.outBoundChoiceTimes);

rewardComponentPatternwise         = zeros(dimension,6);
errorComponentPatternwise          = zeros(dimension,6);
inBoundChoiceComponentPatternwise  = zeros(dimension,6);
outBoundChoiceComponentPatternwise = zeros(dimension,6);

for i = 1:6
    for j = 1:dimension
        try
        rewardComponentPatternwise(j,i) = mean(rewardedComponent(:,(i-1)*5+j));
        catch
            keyboard
        end
        errorComponentPatternwise(j,i) = mean(erroredComponent(:,(i-1)*5+j));
        inBoundChoiceComponentPatternwise(j,i) = mean(inBoundChoiceComponent(:,(i-1)*5+j));
        outBoundChoiceComponentPatternwise(j,i) = mean(outBoundChoiceComponent(:,(i-1)*5+j));
    end
end

X = {rewardComponentPatternwise, errorComponentPatternwise, inBoundChoiceComponentPatternwise, outBoundChoiceComponentPatternwise};

for i = 1:4
    critical_behaviors(i).mean_component_strength = X{i};
end

%%
critical_components = {rewardedComponent, erroredComponent, inBoundChoiceComponent, outBoundChoiceComponent};
critical_times = {rewardedComponentTime, erroredComponentTime, inBoundChoiceComponentTime, outBoundChoiceComponentTime};
for i = 1:4
    critical_behaviors(i).component_strength = critical_components{i};
    critical_behaviors(i).component_time = critical_times{i};
end

%% lindist and tperf, at their own interped times
rewardedLindist       = animal_behavior.lindist(animal_behavior.rewardTimes,:);
erroredLindist        = animal_behavior.lindist(animal_behavior.errorTimes,:);
inBoundChoiceLindist  = animal_behavior.lindist(animal_behavior.inBoundChoiceTimes,:);
outBoundChoiceLindist = animal_behavior.lindist(animal_behavior.outBoundChoiceTimes,:);
L = {rewardedLindist, erroredLindist,inBoundChoiceLindist,outBoundChoiceLindist};

rewardedInboundPerf  = animal_behavior.tperfInbound(animal_behavior.rewardTimes,:);
rewardedOutboundPerf = animal_behavior.tperfOutbound(animal_behavior.rewardTimes,:);

erroredInboundPerf  = animal_behavior.tperfInbound(animal_behavior.errorTimes,:);
erroredOutboundPerf = animal_behavior.tperfOutbound(animal_behavior.errorTimes,:);


inBoundChoicePerf  = animal_behavior.tperfInbound(animal_behavior.inBoundChoiceTimes,:);
outBoundChoicePerf = animal_behavior.tperfOutbound(animal_behavior.outBoundChoiceTimes,:);
P_i = {rewardedInboundPerf,erroredInboundPerf,inBoundChoicePerf ,[]};
P_o = {rewardedOutboundPerf ,erroredOutboundPerf,[],outBoundChoicePerf};


for i = 1:4
    critical_behaviors(i).lindist = L{i};
    critical_behaviors(i).tperfInbound = P_i{i};
    critical_behaviors(i).tperfOutbound = P_o{i};
end

%% tperf, at unified times

% FIND UNIFIED TIMES
tperfInbound = {critical_behaviors.tperfInbound};
tperfOutbound = {critical_behaviors.tperfOutbound};

unifiedTime = [];
if ~ unify_from_running
    for i = 1:numel(critical_times)
        unifiedTime = sort(union(unifiedTime, critical_times{i}));
    end
else
    unifiedTime = animal_behavior.time;
end
% interp
for i = 1:numel(critical_times)
    try
    critical_components{i} = interp1(critical_times{i}, critical_components{i}, unifiedTime, 'linear');
    catch
        keyboard
    end
    if ~isempty(tperfInbound{i})
        tperfInbound{i} = interp1(critical_times{i}, tperfInbound{i}, unifiedTime,'linear');
    end
    if ~isempty(tperfOutbound{i})
        tperfOutbound{i} = interp1(critical_times{i}, tperfOutbound{i}, unifiedTime,'linear');
    end
    critical_times{i} = unifiedTime;
end


for i = 1:4
    critical_behaviors(i).unifiedTperfInbound = tperfInbound{i};
    critical_behaviors(i).unifiedTperfOutbound = tperfOutbound{i};
    
%     [critical_behaviors(i).inBound_learning_times, ...
%         critical_behaviors(i).inBound_stable_times, critical_behaviors(i).inBound_learning_phase, ...
%         critical_behaviors(i).inBound_stable_phase, critical_behaviors(i).inBound_learning_rate, ...
%         critical_behaviors(i).inBound_stable_rate, inbound_throwout_times]...
%         = behaviors.calculateLearning(unifiedTime, tperfInbound{i}, quantile);
%     
%     [critical_behaviors(i).outBound_learning_times, ...
%         critical_behaviors(i).outBound_stable_times, critical_behaviors(i).outBound_learning_phase, ...
%         critical_behaviors(i).outBound_stable_phase, critical_behaviors(i).outBound_learning_rate, ...
%         critical_behaviors(i).outBound_stable_rate, outbound_throwout_times]...
%         =  behaviors.calculateLearning(unifiedTime, tperfOutbound{i}, quantile);
%     
%     [critical_behaviors(i).inBound_stable_component,...
%         critical_behaviors(i).inBound_learning_component, ...
%         critical_behaviors(i).inBound_stable_corr, ...
%         critical_behaviors(i).inBound_learning_corr] = ....
%         components.calculateLearningByComponent(critical_components{i}, inbound_throwout_times,...
%         critical_behaviors(i).inBound_learning_rate, ...
%         critical_behaviors(i).inBound_stable_rate, critical_behaviors(i).inBound_learning_times,...
%         critical_behaviors(i).inBound_stable_times);
%     
%     [critical_behaviors(i).outBound_stable_component,...
%         critical_behaviors(i).outBound_learning_component, ...
%         critical_behaviors(i).outBound_stable_corr, ...
%         critical_behaviors(i).outBound_learning_corr] = ....
%         components.calculateLearningByComponent(critical_components{i}, outbound_throwout_times,...
%         critical_behaviors(i).outBound_learning_rate, ...
%         critical_behaviors(i).outBound_stable_rate, critical_behaviors(i).outBound_learning_times,...
%         critical_behaviors(i).outBound_stable_times);
end
