function t_event = eventuv(s, pattern_labels, event_labels, uv_labels, scalar_info)
    % EVENTSPEC Convert a struct array of event structs to a table

    if nargin < 5 || isempty(scalar_info)
        scalar_info = struct();
    end
    if nargin < 4 || isempty(uv_labels)
        uv_labels = 1:size(s(1).event_u_values, 3);
    end
    if nargin < 3 || isempty(event_labels)
        event_labels = 1:size(s(1).event_u_values, 2);
    end
    if nargin < 2 || isempty(pattern_labels)
        pattern_labels = 1:size(s(1).event_u_values, 1);
    end

    if numel(s) == 1
        % Define the event field names
        event_fields = {'event_u_values', 'event_v_values', 'patterns', 'events', 'uv_components'};

        s.patterns = pattern_labels(:);
        s.events = event_labels(:)';
        s.uv_components = permute(uv_labels(:), [3, 2, 1]);
        s = rmfield(s, setdiff(fieldnames(s), event_fields));
        s = nd.broadcast(s);

        % Initialize the cell arrays for the event columns
        event_columns = cell(1, numel(event_fields));  
        event_names = event_fields;

        % Add the event columns
        for i = 1:numel(event_fields)
            field = event_fields{i};
            array = reshape(s.(field), [], 1);
            event_columns{i} = array;
        end

        % Add scalar info to the table
        scalar_fields = fieldnames(scalar_info);
        for i = 1:numel(scalar_fields)
            field = scalar_fields{i};
            scalar_value = scalar_info.(field);
            event_columns{end+1} = repmat(scalar_value, [numel(event_columns{1}), 1]);
            event_names{end+1} = field;
        end
        
        % Check that all columns have the same number of rows
        event_lengths = cellfun(@numel, event_columns);
        
        if range(event_lengths) ~= 0
            error('The event columns do not all have the same number of rows.')
        end

        % Construct the table
        t_event = table(event_columns{:}, 'VariableNames', event_names);
    else
        t_event = cell(1, numel(s));
        for i = 1:numel(s)
            if nd.isEmpty(s(i))
                continue
            end
            t_event{i} = table.analyses.eventuv(s(i), pattern_labels, event_labels, uv_labels, scalar_info);
        end
        t_event(cellfun(@isempty, t_event)) = [];
        t_event = vertcat(t_event{:});
    end
end

