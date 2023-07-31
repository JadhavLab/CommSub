function [t_uv, t_spec] = trigspec(s, time_labels, frequency_labels, scalar_info)
    % TRIGSPEC Convert a struct array of trigspec structs to two tables
    if numel(s) > 1
        t_uv = cell(1, numel(s));
        t_spec = cell(1, numel(s));
        for i = 1:numel(s)
            [t_uv{i}, t_spec{i}] = table.analyses.trigspec(s(i), time_labels, frequency_labels, scalar_info);
        end
        t_uv   = vertcat(t_uv{:});
        t_spec = vertcat(t_spec{:});
    else

        s = nd.unnest(s, 'spec_avg');
        s = nd.broadcast(s);
        s = rmfield(s, ["threshold_crossed_times", "spec_stderr"]);
        if nargin < 4
            scalar_info = struct();
        end
        if nargin < 3
            frequency_labels = 1:size(s(1).S1, 2);
        end
        if nargin < 2
            time_labels = 1:size(s(1).S1, 1);
        end

        % Define the UV and spectral component field names
        uv_fields = {'u_average', 'v_average', 'u_stderr', 'v_stderr'};
        common_fields = {'name', 'comp', 'direction'};
        spec_fields = setdiff(fieldnames(s), [uv_fields, common_fields]);

        % Initialize the cell arrays for the UV and spectral component columns
        uv_columns = cell(1, numel(uv_fields) + numel(common_fields) + 2);
        uv_names = [uv_fields(:)', common_fields(:)', {'time', 'components'}];
        
        spec_columns = cell(1, numel(spec_fields) + numel(common_fields) + 2);
        spec_names = [spec_fields(:)', common_fields(:)', {'time', 'frequency'}];

        % Add the UV component columns
        for i = 1:numel(uv_fields)
            field = uv_fields{i};
            array = reshape(s.(field), [], 1);
            uv_columns{i} = array;
        end

        % Add the spectral component columns
        for i = 1:numel(spec_fields)
            field = spec_fields{i};
            array = reshape(s.(field), [], 1);
            spec_columns{i} = array;
        end

        % Add the common fields to both tables
        for i = 1:numel(common_fields)
            field = common_fields{i};
            array_uv = reshape(s.(field)(:, 1:size(s.u_average, 2)), [], 1);
            array_spec = reshape(s.(field), [], 1);
            uv_columns{numel(uv_fields) + i} = array_uv;
            spec_columns{numel(spec_fields) + i} = array_spec;
        end

        % Create the time and frequency labels
        uv_columns{end-1} = repmat(time_labels(:), [size(s.u_average, 2), 1]);
        uv_columns{end}   = reshape(repmat((1:size(s.u_average, 2)), [numel(time_labels), 1]), [], 1);

        spec_columns{end-1} = repmat(time_labels(:), [size(s.(spec_fields{1}), 2), 1]);
        spec_columns{end} = reshape(repmat(frequency_labels(:)', [numel(time_labels), 1]), [], 1);

        % Add scalar info to each table
        scalar_fields = fieldnames(scalar_info);
        for i = 1:numel(scalar_fields)
            field = scalar_fields{i};
            scalar_value = scalar_info.(field);
            uv_columns{end+1} = repmat(scalar_value, [numel(uv_columns{1}), 1]);
            uv_names{end+1} = field;
            spec_columns{end+1} = repmat(scalar_value, [numel(spec_columns{1}), 1]);
            spec_names{end+1} = field;
        end
        
        % Check that all columns have the same number of rows
        uv_lengths   = cellfun(@numel, uv_columns);
        spec_lengths = cellfun(@numel, spec_columns);
        
        if range(uv_lengths) ~= 0
            error('The UV columns do not all have the same number of rows.')
        end
        
        if range(spec_lengths) ~= 0
            error('The spectral columns do not all have the same number of rows.')
        end

        % Construct the tables
        t_uv   = table(uv_columns{:}, 'VariableNames', uv_names);
        t_spec = table(spec_columns{:}, 'VariableNames', spec_names);

    end

end

