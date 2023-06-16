function [trajRewardTimes, trajErrorTimes, trajChoiceTimes] = findTrajCriticalTimes(animal)
% Finds the rewarded and error times from trajinfo file
% Finds the choice times from the lregion file
% where during outbound 0.2-0.4 are used and during inbound, 0.8-0.7 and
% 0.6-0.5 are used

% Output:
% trajRewardTimes: cell that contains the inbound(1) and outbound(2)
% rewarded times
% trajErrorTimes: cell that contains the inbound(1) and outbound(2)
% errored times
% trajChoiceTimes: cell that contains the inbound(1) and outbound(2)
% decision times
% otherTimes: cell that contains less critical times during inbound(1) and outbound(2)
%% rewarded and errored times
load(animal + "trajinfo01.mat");
trajinfo = ndb.toNd(trajinfo);
trajRewardTimes = cell(1,2);
trajErrorTimes = cell(1,2);
trajChoiceTimes = cell(1,2);

trajBoundTimes = cell(1,2);
otherTimes = cell(1,2);

for i = 2:2:numel(trajinfo)
    rewarded = trajinfo(i).rewarded;
    trajBound = trajinfo(i).trajbound;

    trajRewardTimes{1} = [trajRewardTimes{1}, trajinfo(i).trajtime(rewarded==1&trajBound==1,2)'];
    trajRewardTimes{2} = [trajRewardTimes{2}, trajinfo(i).trajtime(rewarded==1&trajBound==0,2)'];
    trajErrorTimes{1} = [trajErrorTimes{1}, trajinfo(i).trajtime(rewarded==1&trajBound==1,2)'];
    trajErrorTimes{2} = [trajErrorTimes{2}, trajinfo(i).trajtime(rewarded==0&trajBound==0,2)'];
end

%% choice times
% TR [trajectory, track_region]
% TB [rewarded, trajbound]

load(animal + "lregion01.mat");
lregion = ndb.toNd(lregion);
inBoundChoiceTimes = [];
outBoundChoiceTimes = [];

% find the inbound and outbound choice times
for i = 2:2:numel(lregion)
    trajBound = lregion(i).TB(:,2);
    times = lregion(i).time;
    trajBoundTimes{1} = [trajBoundTimes{1}, times(trajBound == 1)'];
    trajBoundTimes{2} = [trajBoundTimes{2}, times(trajBound == 0)'];
    
    linDist = lregion(i).lindist;
    inBoundChoiceTimes = [inBoundChoiceTimes, ...
                          times(trajBound == 1 & linDist >= 0.2 & linDist <= 0.4)'];
    outBoundChoiceTimes = [outBoundChoiceTimes, times(trajBound == 0 & ...
                          ((linDist >= 0.7 & linDist <= 0.8)| (linDist >= 0.5 & linDist <= 0.6)))'];
end

trajChoiceTimes{1} = inBoundChoiceTimes;
trajChoiceTimes{2} = outBoundChoiceTimes;


otherTimes{1} = setdiff(trajBoundTimes{1}, trajChoiceTimes{1});
otherTimes{1} = setdiff(otherTimes{1}, trajRewardTimes{1});
otherTimes{1} = setdiff(otherTimes{1}, trajErrorTimes{1});


otherTimes{2} = setdiff(trajBoundTimes{2}, trajChoiceTimes{2});
otherTimes{2} = setdiff(otherTimes{2}, trajRewardTimes{2});
otherTimes{2} = setdiff(otherTimes{2}, trajErrorTimes{2});
end

