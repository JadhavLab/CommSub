function savedat(Patterns_overall, behavior, Option, figAppend)
    [U, V] = deal(Patterns_overall(end).cca.u, Patterns_overall(end).cca.v);
    uv_time = Patterns_overall(end).X_time;
    % Convert behavior table to pandas DataFrame
    behavior_dict = table2struct(behavior, 'ToScalar', true);
    save(figuredefine("data", figAppend + ".mat"), 'U', 'V', 'uv_time', 'Option', 'behavior_dict');
end

