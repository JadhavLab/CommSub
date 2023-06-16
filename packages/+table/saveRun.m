function [runGroupSummaries, individualPartitions, hash] = ...
                                    saveRun(Patterns, Raw, Option)
% Packages data into tables and determinees unique hash identity for a run
% with this set of options. If a run has been previously executed with the has
% the table entries overwritee the previous matching entries.

% For Ziyi's station
path = "C:\Users\BrainMaker\MATLAB Drive\Shared";
addpath(path)

% Prep our two tables
% -------------------
if exist(table.filename.individualPartitions(), 'file')
    load(table.filename.individualPartitions());
else
    runGroupSummaries = [];
end
if exist(table.filename.runGroupSummaries(), 'file')
    load(table.filename.runGroupSummaries());
else
    warning("no exisiting real table")
    individualPartitions = [];
end

% Determine how optimal our run was
% ---------------------------------
[results] = params.optimizeOptionsScript(Patterns, Raw, Options);

Optiontable  = struct2table(Option, 'AsArray', true);
Patterntable = query.getPatternTable(Patterns);

if sum(Patterntable.singularWarning) > 0.7 * height(Patterntable)
    warning("too many singular results in factor analysis")
    choice = input ("proceed to store results?");
    if choice == "no"
        error ("unsucessful run")
    end
end

% Identifying information about this options set and date of run
hash = DataHash(Option);
hash = hash(1:7); % Take the first 7 letters of the hash
hash = string(hash);
timestamp = string(date());

nWindowsCut = windows.getNumWindowsCut(cellOfWindows);
optimizationResult = result.optimizationResult;
summaryAddon = table(timestamp, hash, getNumWindowsCut, cutoffs, ...
                     );
thisRunSummary = [Optiontable, summaryAddon]; % Combine option columnns with hash and date

partitionAddon = repmat(thisRunSummary, height(Patterntable), 1);
thisRunsPartitions = [Patterntable, partitionAddon]; % combine those with all rows of the pattern table

% Check and Hash
if ~isempty('runGroupSummaries') && any(contains(TABLE.hash, hash))

    matchingHashLocation = contains(runGroupSummaries.hash, hash);
    runGroupSummaries(matchingHashLocation, :) = []; 

    matchingHashLocation = contains(individualPartitions.hash, hash);
    individualPartitions(matchingHashLocation, :) = [];

    disp("already computed before, rehashing to the same location");
else

    disp("new results stored!")

end
old_height = height(individualPartitions);


% Append the new results
newColumnsWereAdded = width(individualPartitions) ~= width(thisRunsPartitions);
if newColumnsWereAdded
    individualPartitions = table.addNewColumn(individualPartitions, thisRunsPartitions);
    runGroupSummaries    = table.addNewColumn(runGroupSummaries,    thisRunSummary);
else
    individualPartitions = [individualPartitions; thisRunsPartitions];
    runGroupSummaries    = [runGroupSummaries;    thisRunSummary];
end
assert(height(individualPartitions) > old_height, "not appending");

% Save tables
save(table.filename.runGroupSummaries(),    "runGroupSummaries",    '-v7.3');
save(table.filename.individualPartitions(), "individualPartitions", '-v7.3');

