function [animal_behavior, unique_times] = addBehToTable(animal_behavior)
% [animal_behavior, unique_times] = addBehToTable(animal_behavior)
%
% ADDBEHTOTABLE add some more entries to the animal behavioral table
%
% input: 
% - animal_behavior: the animal behavioral table
%
% output: 
% - animal_behavior: the animal behavioral table with more entries
% - unique_times: the unique times in the table

[~, unique_times, ~] = unique(animal_behavior.time);
animal_behavior = animal_behavior( unique_times,:);

%% identify the decision, error, reward times to the table
animal_behavior.wellTimes = (animal_behavior.trajbound == 0  & animal_behavior.lindist > 0.99)...
                            |(animal_behavior.trajbound == 1 & animal_behavior.lindist <0.01) ; %0.01 inbound
animal_behavior.rewardTimes = animal_behavior.wellTimes & animal_behavior.rewarded ==1;
animal_behavior.outBounderrorTimes = animal_behavior.wellTimes & animal_behavior.rewarded ==0;
animal_behavior.inBounderrorTimes = animal_behavior.trajbound == 1 & animal_behavior.lindist > 0.99  & animal_behavior.rewarded ==0;
animal_behavior.errorTimes = animal_behavior.outBounderrorTimes | animal_behavior.inBounderrorTimes;

animal_behavior.outBoundChoiceTimes = animal_behavior.trajbound == 0 & animal_behavior.lindist >= 0.2 & animal_behavior.lindist <= 0.4;
animal_behavior.inBoundChoiceTimes = animal_behavior.trajbound == 1 & animal_behavior.lindist >= 0.5 & animal_behavior.lindist <= 0.6;

%% interp the inbound and outbound performances to all times
tperfInbound = (animal_behavior.tperf_timewise(animal_behavior.trajbound ==1));
tperfOutbound = (animal_behavior.tperf_timewise(animal_behavior.trajbound ==0));

animal_behavior.tperfInbound = interp1(animal_behavior.time(animal_behavior.trajbound ==1),...
               tperfInbound, animal_behavior.time);
           
animal_behavior.tperfOutbound = interp1(animal_behavior.time(animal_behavior.trajbound==0),...
               tperfOutbound, animal_behavior.time);

end

