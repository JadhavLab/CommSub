function [subspaceDist, rowVar] = sequentialDistanceMatrix(rt, K)

if nargin < 2
    K = inf;
end

rt = rt(rt.K <= K, :);

% Get unique combinations of baseDirection, basePattern, and method
groups = unique(rt(:, ["baseDirection", "basePattern", "method"]));
groups2 = unique(rt(:, ["removeDirection", "removePattern", "dimensionRemoved"]));

% Initialize a cell array to hold the performance vectors
performanceVectors = cell(height(groups), 1);

% Loop over unique combinations
for i = 1:height(groups)
    % Extract the rows of rt that match the current combination
    mask = ismember(rt(:, ["baseDirection", "basePattern", "method"]), groups(i, :));
    subrt = rt(mask, :);

    % Create a vector of performance values for each unique combination of removeDirection, removePattern, and dimensionRemoved
    subgroups = groups2;
    performanceVector = nan(height(subgroups), 1);
    for j = 1:height(subgroups)
        mask = ismember(subrt(:, ["removeDirection", "removePattern", "dimensionRemoved"]), subgroups(j, :));
        performance = subrt.performance(mask);
        if isempty(performance)
            performanceVector(j) = 0;
        else
            performanceVector(j) = performance(1);
        end
    end

    % Store the performance vector
    performanceVectors{i} = performanceVector;
end

cellfun(@numel, performanceVectors)

% Calculate the distance matrix
subspaceDist = nan(height(groups), height(groups));
for i = 1:height(groups)
    for j = 1:height(groups)
        if i == j
            continue
        end
        subspaceDist(i, j) = pdist2(performanceVectors{i}(:)', performanceVectors{j}(:)', "euclidean")
    end
end

rowVar = string(nan(height(groups), 1));
directs = ["hH", "hP"];
for i = 1:height(groups)
    rowVar(i) = directs(groups(i, :).baseDirection) + "_" + groups(i, :).basePattern + "_"+  groups(i, :).method;
end
