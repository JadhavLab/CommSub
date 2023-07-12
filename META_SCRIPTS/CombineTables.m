% This script combines all the runs summary files into one table

%% Combine all the RunSummary files into one table
regex = fullfile(codedefine, "DATA_TABLES", "RunsSummary_*.mat");
finalFile = "RunsSummary.mat";
out = table.combineAndUpdateTables(regex, finalFile);
disp(head(out,10))

%% Comine the DetailedRunsSummary files into one table
regex = fullfile(codedefine, "DATA_TABLES", "DetailedRunsSummary_*.mat");
finalFile = "DetailedRunsSummary.mat";
out = table.combineAndUpdateTables(regex, finalFile);
disp(head(out,10))

