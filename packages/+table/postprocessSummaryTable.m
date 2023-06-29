function [summaryTable] = postprocessSummaryTable(summaryTable)
% POSTPROCESSSUMMARYTABLES Postprocess summary tables
%
% This first reorders the columns, so that timestamp and hash are first, and
% then sorts by timestamp.  We keep do not toss out hash and timestamp.columns

% Column names as strings
colNames = string(summaryTable.Properties.VariableNames);
% Setdiff to get the columns that are not timestamp or hash
nonTimestampOrHash = setdiff(colNames, ["timestamp", "hash"]);
% New order
newOrder = ["timestamp", "hash", nonTimestampOrHash];
% Reorder
summaryTable = summaryTable(:, newOrder);

% Sort by timestamp
summaryTable = sortrows(summaryTable, "timestamp");

summaryTable = util.table.castefficient(summaryTable);

end
