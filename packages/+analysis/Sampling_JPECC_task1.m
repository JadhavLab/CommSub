function [tens_curr_source, tens_curr_target] = Sampling_JPECC_task1(Source, Target)

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
    selected_source_idx = zeros(1, num_samples_target);

    % For each of the first n target times
    for d = 1:num_samples_target
        % Calculate the time differences between the current target time and all source times
        time_diffs = abs(source_time - target_time(d));
    
        % Find the source time that is closest to the current target time
        [~, closest_source_idx] = min(time_diffs);
    
        % Add the index of the closest source time to the selected source indices
        selected_source_idx(d) = source_idx(closest_source_idx);
    
        % Remove the closest source time from consideration in future iterations
        source_time(closest_source_idx) = [];
        source_idx(closest_source_idx) = [];
    end

    tens_curr_source = tens_curr_source(selected_source_idx, :);

elseif num_samples_source < num_samples_target

    % Sort the time vectors
    [source_time, ~] = sort(source_time);
    [target_time, target_idx] = sort(target_time);

    % Initialize the selected source indices
    selected_target_idx = zeros(1, num_samples_source);

    % For each of the first n target times
    for d = 1:num_samples_source
        % Calculate the time differences between the current target time and all source times
        time_diffs = abs(target_time - source_time(d));
    
        % Find the source time that is closest to the current target time
        [~, closest_target_idx] = min(time_diffs);
    
        % Add the index of the closest source time to the selected source indices
        selected_target_idx(d) = target_idx(closest_target_idx);
    
        % Remove the closest source time from consideration in future iterations
        target_time(closest_target_idx) = [];
        target_idx(closest_target_idx) = [];
   end

   tens_curr_target = tens_curr_target(selected_target_idx, :);

end