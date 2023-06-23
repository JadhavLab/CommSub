
source = Patterns.X_source;
target = Patterns.X_target;
% source_index = Patterns.index_source;
% target_index = Patterns.index_target;
both = [source;target];

rr = Patterns.rankRegress
B_ = rr.B_;
[u,~,v] = svd(B_);

dips("Size source: " + size(source))
disp("Size u: " + size(u))
disp("Size target: " + size(target))
disp("Size v: " + size(v))


match_inputoutput = [u(:,1:3);v(:,1:3)];
size(match_inputoutput), size(both)

activities = match_inputoutput' * both;

figure
plot(abs(activities'))
legend({'component 1', 'component 2', 'component 3'})
title('Matched input-output components')

% figure
plots.event.plotSingleEvents(cellOfWindows, Events, Option, r, ...
    'before', 0.5, 'after', 0.5, 'appendFigTitle', char(Option.animal), ...
    'eventTypes', ["theta", "delta", "ripple"], ...
    'subspaceComponents', activities, ...
    'savePath', fullfile(figuredefine, 'singleEvent_plusActivities'), ...
    'displayType', 'heatmap');
