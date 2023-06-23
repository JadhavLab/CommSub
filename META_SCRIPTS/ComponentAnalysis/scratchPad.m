
source = Patterns.X_source;
target = Patterns.X_target;
time   = Patterns.X_time;
source_index = Patterns.index_source;
target_index = Patterns.index_target;
source = r.spikeCountMatrix(source_index,:);
target = r.spikeCountMatrix(target_index,:);
time = r.timeBinMidPoints;
both = [source;target];

rr = Patterns.rankRegress
B_ = rr.B_;
[u,s,v] = svd(B_, 'full');
sdiag = diag(s);

dips("Size source: " + size(source))
disp("Size u: " + size(u))
disp("Size target: " + size(target))
disp("Size v: " + size(v))

% s_u = s(1:3,1:3);
% s_v = s(1:3,1:3);
s_u = 1;
s_v = 1;
match_inputoutput = [u(:,1:3)*s_u;v(:,1:3)*s_v];
size(match_inputoutput), size(both)
activities = match_inputoutput' * (both);
activities = activities - mean(activities,2);
activities = abs(activities);

figure
plot(abs(activities'))
legend({'component 1', 'component 2', 'component 3'})
title('Matched input-output components')

% figure
plots.event.plotSingleEvents(cellOfWindows, Events, Option, r, ...
    'before', 0.5, 'after', 0.5, 'appendFigTitle', char(Option.animal), ...
    'eventTypes', ["theta", "delta", "ripple"], ...
    'subspaceComponents', activities, ...
    'subspaceTimes', time, ...
    'saveFolder', fullfile(figuredefine, 'singleEvent_plusActivities'), ...
    'displayType', 'heatmap');
