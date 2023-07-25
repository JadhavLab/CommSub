function new_table = flexvertcat(varargin)
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


if numel(varargin) == 2
    A = varargin{1};
    B = varargin{2};
    %% find the differences in column names
    varsA = A.Properties.VariableNames;
    varsB = B.Properties.VariableNames;

    %% Check if the columns that exist in both tables have the same type.
    common_columns = intersect(varsA, varsB);
    for i = 1:numel(common_columns)
        columnName = common_columns{i};
        if ~isa(A.(columnName), class(B.(columnName)))
            try
                A.(columnName) = feval(class(B.(columnName)), A.(columnName));
            catch ME
                if isa(A.(columnName), 'cell') && all(cellfun(@isempty, A.(columnName)))
                    A.(columnName) = repmat(feval(class(B.(columnName)), NaN), height(A), 1);
                else
                    warning('Failed to convert column %s from %s to %s: %s', ...
                            columnName, class(A.(columnName)), class(B.(columnName)), ME.message);
                end
            end
        end
    end

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
    try
        [A, B] = table.alignTables(A, B, 'aggressive', true);
        new_table = [A; B];
    catch ME1
        % try
            new_table = aggressive_match(A, B);
        % catch ME2
        %     keyboard
        % end
    end
    return;

else
    new_table = varargin{1};
    for i = 2:numel(varargin)
        new_table = table.flexvertcat(new_table, varargin{i});
    end
end

function C = aggressive_match(A, B)
    % Assuming that 'A' and 'B' are your tables
    % Get sizes of A and B
    [~, num_cols_A] = size(A);
    [~, num_cols_B] = size(B);

    if num_cols_A ~= num_cols_B
        error('The tables have different number of columns and cannot be concatenated.')
    else
        % Iterate through each column and check the data type
        F = setdiff(fieldnames(A), {'Properties', 'Row', 'VarNames', 'Variables'});
        F = F(:)';
        for f = F
            col_A = A.(f{1});
            col_B = B.(f{1});
            try
                [col_A; col_B];
                concatable = true;
            catch
                concatable = false;
            end
            
            % If the data types do not match
            if ~concatable || ~strcmp(class(col_A), class(col_B))
                of1 = A.Properties.VariableNames;
                of2 = B.Properties.VariableNames;
                assert(numel(setdiff(of1, of2))==0, ...
                    'The tables have different column names and cannot be concatenated.')
                % Remove the column from the table
                A(:, f{1}) = [];
                B(:, f{1}) = [];
                % Convert the columns to strings in both tables
                try
                    col_A = string(col_A);
                    col_B = string(col_B);
                catch
                    continue
                end
                ncolsA = size(col_A(1,:));
                ncolsB = size(col_B(1,:));
                if ~isequal(ncolsA, ncolsB)
                    if any(ncolsA > ncolsB)
                        col_B = repmat(col_B, 1, ncolsA(2) - ncolsB(2) + 1);
                    elseif any(ncolsA < ncolsB)
                        col_A = repmat(col_A, 1, ncolsB(2) - ncolsA(2) + 1);
                    end
                end
                A.(f{1}) = col_A;
                B.(f{1}) = col_B;
                % Set order
                A = A(:, of1);
                B = B(:, of2);
            end
        end
    end
    % Concatenate the tables
    C = [A; B];
