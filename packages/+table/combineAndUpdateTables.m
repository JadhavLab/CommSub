function out = combineAndUpdateTables(regex, finalFile)
% COMBINEANDUPDATETABLES Combine and update tables
%   COMBINEANDUPDATETABLES(regex, finalFile) combines and updates tables
%   that match the regular expression regex and saves the combined table to
%   finalFile. The tables must have a hash column.
%
%   Example:
%       combineAndUpdateTables('*.mat', 'finalTable.mat')

    % Load the existing table if it exists
    if exist(finalFile, 'file')
        load(finalFile);
    else
        out = [];
    end

    % Get a list of the new table files
    newFiles = dir(regex);
    
    % Loop over each new table file
    for i = 1:length(newFiles)

        % Load the new table
        newTable = load(fullfile(newFiles(i).folder));
        fields = fieldnames(newTable);
        assert(numel(fields) == 1);
        newTable = newTable.(fields{1});

        % Align the existing table and the new table
        [out, newTable] = table.alignTables(out, newTable);
        
        % Loop over each row in the new table
        for j = 1:height(newTable)
            % Find if the hash value exists in the existing table
            idx = find(strcmp(out.hash, newTable.hash(j)));
            
            if isempty(idx)
                % If the hash value does not exist, append the new row
                out = [out; newTable(j, :)];
            else
                % If the hash value exists, replace the existing row
                out(idx, :) = newTable(j, :);
            end
        end
    end
    
    % Save the updated table
    save(finalFile, 'out');
end
