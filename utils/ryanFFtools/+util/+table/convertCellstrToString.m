function tbl = convertCellstrToString(tbl)
%CONVERTCELLSTRTOSTRING Convert cell array of character vectors to string array

    % Get variable names
    varNames = tbl.Properties.VariableNames;

    % Loop over each variable
    for i = 1:length(varNames)
        var = tbl.(varNames{i});
        % Check if variable is a cell array of character vectors
        if iscell(var) && ischar(var{1})
            % Convert to string array
            tbl.(varNames{i}) = string(var);
        end
    end
end

