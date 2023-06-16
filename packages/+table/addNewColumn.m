function new_table = addNewColumn(old_table, new_subTable)
% this function adds a subtable with one or more new columns to the
% original table. The newly-added fields will be set as NaN for the entries
% in the old table

% input:
% old_table: the original table
% new_subTabel: the table to be appended, having more columns

% output:
% new_table


%% find the differences in column names
old_columns = old_table.Properties.VariableNames;
new_columns = new_subTable.Properties.VariableNames;

if numel(old_columns) < numel(new_columns)
    add_columns = setdiff(new_columns, old_columns);
    new_table = old_table;
    for i = 1:numel(add_columns)
        new_columnName = add_columns{i};
        new_table.(new_columnName) = nan(height(new_table),1);
    end
    new_table = [new_table; new_subTable];
else
    add_columns = setdiff(old_columns, new_columns);
    new_table = new_subTable;
    for i = 1:numel(add_columns)
        new_columnName = add_columns{i};
        new_table.(new_columnName) = nan(height(new_table),1);
    end
    new_table = [old_table; new_table];
    
end



end

