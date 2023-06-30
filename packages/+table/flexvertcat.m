function new_table = flexvertcat(A, B)
%ADDNEWCOLUMN add a new column to a table
%
% this function adds a subtable with one or more new columns to the original
% table. The newly-added fields will be set as NaN for the entries in the old
% table or vice versa.
%
% input:
% A: the original table
% B: the subtable to be added
%
% output:
% new_table


%% find the differences in column names
varsA = A.Properties.VariableNames;
varsB = B.Properties.VariableNames;

%% determine columns in A but not in B
add_columns = setdiff(varsA, varsB);
for i = 1:numel(add_columns)
    new_columnName = add_columns{i};
    n = numel(A.(new_columnName)(1,:));
    B.(new_columnName) = nan(height(B), n);
end
%% determine columns in B but not in A
add_columns = setdiff(varsB, varsA);
for i = 1:numel(add_columns)
    new_columnName = add_columns{i};
    n = numel(B.(new_columnName)(1,:));
    A.(new_columnName) = nan(height(A),n);
end
%% Concatenate the two tables
new_table = [A; B];
return;

% if numel(old_columns) < numel(new_columns)
%     add_columns = setdiff(new_columns, old_columns);
%     new_table = old_table;
%     for i = 1:numel(add_columns)
%         new_columnName = add_columns{i};
%         new_table.(new_columnName) = nan(height(new_table),1);
%     end
%     % new_table = [new_table; new_subTable];
%     new_table = util.cell.icat({new_table, new_subTable}, ...
%                                 'fieldCombine', 'union');
% else
%     add_columns = setdiff(old_columns, new_columns);
%     new_table = new_subTable;
%     for i = 1:numel(add_columns)
%         new_columnName = add_columns{i};
%         new_table.(new_columnName) = nan(height(new_table),1);
%     end
%     % new_table = [old_table; new_table];
%     new_table = util.table.icat({old_table, new_subTable}, ...
%                                 'fieldCombine', 'union');
%     
% end



end

