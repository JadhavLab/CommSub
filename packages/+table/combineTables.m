function [newTable, tables] = combineTables(regex, newTableName)
% combineTables - combines all the tables matching a regex into one table
% with a given name
%
% Parameters:
%   regex - the regex to match the tables to combine
%   newTableName - the name of the new table
%
% Returns:
%   newTable - the new table, also saved to a file at the regexp location

% get all the files matching the regex
files = dir(regex);
if isempty(files)
    error('No files matching the regex');
end
% get the number of files
numFiles = length(files);
% create a cell array to hold the tables
tables = cell(numFiles, 1);
% loop through the files
for i = 1:numFiles
    % load the matfiles with the tables
    tmp = load(files(i).name);
    f = fieldnames(tmp);
    assert(numel(f) == 1, 'Too many fields in the matfile');
    tables{i} = tmp.(f{1});
end
% combine the tables, two at a time, with table.flexvertcat
newTable = tables{1};
for i = 2:numFiles
    newTable = table.flexvertcat(newTable, tables{i});
end
% save the table
out.(replace(newTableName, '.mat', '')) = newTable;
if contains(newTableName, filesep)
    save(newTableName, '-struct', 'out', '-v7.3');
else
    [folder,~] = fileparts(regex);
    save(char(fullfile(folder, newTableName)), '-struct', 'out', '-v7.3');
end
