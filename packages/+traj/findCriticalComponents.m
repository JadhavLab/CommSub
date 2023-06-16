function [rewardComponents, errorComponents, choiceComponents, otherComponents] = ...
                             findCriticalComponents(running_subspaces, ...
                spikeTime,trajRewardTimes, trajErrorTimes, trajChoiceTimes, otherTimes)

% this function returns the components corresponding to the critical times,
% as well as the components during less critical times along the track

% output

%% find the corresponding times in the subspace matrix
rewardComponents = cell(1,2);
errorComponents = cell(1,2); 
choiceComponents = cell(1,2); 
otherComponents = cell(1,2);


for i = 1:2
indicesTimesRewardComponents = BetweenTimes(spikeTimes,trajRewardTimes{i});
indicesTimesErrorComponents = BetweenTimes(spikeTimes,trajErrorTimes{i});
indicesTimesChoiceComponents = BetweenTimes(spikeTimes,trajChoiceTimes{i});
indicesTimesOtherComponents = BetweenTimes(spikeTimes,otherTimes{i});


rewardComponents{i} = running_subspaces(:,indicesTimesRewardComponents(~isnan(indicesTimesRewardComponents)));
errorComponents{i} = running_subspaces(:,indicesTimesErrorComponents(~isnan(indicesTimesErrorComponents)));
choiceComponents{i} = running_subspaces(:,indicesTimesChoiceComponents(~isnan(indicesTimesChoiceComponents)));
otherComponents{i} = running_subspaces(:,indicesTimesOtherComponents(~isnan(indicesTimesOtherComponents)));

end

