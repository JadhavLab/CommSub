function [tens_curr_source, tens_curr_target, is_source] = Sampling_JPECC_task2(Source, Target, upper, lower)
%is_source = 1/0, represent whether tens_curr_target is a 1*3 cell

tens_curr_source = Source.val;
tens_curr_target = Target.val;
source_time = Source.time;
target_time = Target.time;

% sampling
num_samples_source = size(tens_curr_source, 1);
num_samples_target = size(tens_curr_target, 1);

if num_samples_source > num_samples_target

    % Sort the time vectors
    [source_time, source_idx] = sort(source_time);
    [target_time, ~] = sort(target_time);

    % Initialize the selected source indices
    selected_source_idx = cell(1, num_samples_target);

    % For each of the first n target times
    for d = 1:num_samples_target
        % Calculate the time differences between the current target time and all source times
        time_diffs = source_time - target_time(d);
        
        % Find the source times that are within the threshold to the current target time
        close_source_idx = find(time_diffs <= upper & time_diffs >= lower);
        
        % If there are more than 5 close source times, randomly select 5 of them
        if length(close_source_idx) > 5
            close_source_idx = randsample(close_source_idx, 5);
        end
        
        % Add the indices of the close source times to the selected source indices
        selected_source_idx{d} = source_idx(close_source_idx);
        
        % Remove the close source times from consideration in future iterations
        source_time(close_source_idx) = [];
        source_idx(close_source_idx) = [];
    end

    % Remove empty cells in selected_source_idx and corresponding elements in target
    empty_cells = cellfun('isempty', selected_source_idx);
    selected_source_idx(empty_cells) = [];
    tens_curr_target = tens_curr_target(~empty_cells,:);

    % Initialize a cell array to store the three sets of tens_curr_source
    tens_curr_source_sets = cell(1, 30);
    
    % Randomly select one source time from each cell in selected_source_idx three times
    for i = 1:30
        selected_source_idx_sample = cellfun(@(x) x(randi(length(x))), selected_source_idx);
        
        % Check for duplicates in selected_source_idx_sample and resample if necessary
        while length(unique(selected_source_idx_sample)) < length(selected_source_idx_sample)
            duplicate_idx = find(hist(selected_source_idx_sample, unique(selected_source_idx_sample)) > 1);
            for j = 1:length(duplicate_idx)
                selected_source_idx_sample(duplicate_idx(j)) = randsample(selected_source_idx{duplicate_idx(j)}, 1);
            end
        end
        
        tens_curr_source_sets{i} = tens_curr_source(selected_source_idx_sample, :);
    end
        tens_curr_source = tens_curr_source_sets;
        is_source = 1;

elseif num_samples_source < num_samples_target

    % Sort the time vectors
    [source_time, ~] = sort(source_time);
    [target_time, target_idx] = sort(target_time);

    % Initialize the selected target indices
    selected_target_idx = cell(1, num_samples_source);

    % For each of the first n source times
    for d = 1:num_samples_source
        % Calculate the time differences between the current source time and all target times
        time_diffs = source_time(d) - target_time;
        
        % Find the target times that are within the threshold to the current source time
        close_target_idx = find(time_diffs <= upper & time_diffs >= lower);
        
        % If there are more than 5 close target times, randomly select 5 of them
        if length(close_target_idx) > 5
            close_target_idx = randsample(close_target_idx, 5);
        end
        
        % Add the indices of the close target times to the selected target indices
        selected_target_idx{d} = target_idx(close_target_idx);
        
        % Remove the close target times from consideration in future iterations
        target_time(close_target_idx) = [];
        target_idx(close_target_idx) = [];
    end

    % Remove empty cells in selected_target_idx and corresponding elements in source
    empty_cells = cellfun('isempty', selected_target_idx);
    selected_target_idx(empty_cells) = [];
    tens_curr_source = tens_curr_source(~empty_cells,:);
    
    % Initialize a cell array to store the three sets of tens_curr_source
    tens_curr_target_sets = cell(1, 30);
    
    % Randomly select one source time from each cell in selected_source_idx three times
    for i = 1:30
        selected_target_idx_sample = cellfun(@(x) x(randi(length(x))), selected_target_idx);
        
        % Check for duplicates in selected_source_idx_sample and resample if necessary
        while length(unique(selected_target_idx_sample)) < length(selected_target_idx_sample)
            duplicate_idx = find(hist(selected_target_idx_sample, unique(selected_target_idx_sample)) > 1);
            for j = 1:length(duplicate_idx)
                selected_target_idx_sample(duplicate_idx(j)) = randsample(selected_target_idx{duplicate_idx(j)}, 1);
            end
        end
        
        tens_curr_target_sets{i} = tens_curr_target(selected_target_idx_sample, :);
    end
        tens_curr_target = tens_curr_target_sets;
        is_source = 0;
else
    is_source = 2;

end
