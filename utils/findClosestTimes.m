function [indices, times_found] = findClosestTimes(original_series, search_series)

% given two series of time, find the indices of times in the search series
% that most approximates the times in the original series

tracker = 1; % tracks where in the search series the last found time point 
             % was, so that less reptitive search would be performed
indices = zeros(1,numel(original_series));
is_indices = zeros(1,numel(original_series));
%%
% search the larger series for the closest time points and return their
% indices within the larger array
for i = 1:numel(indices)
    tmp = original_series(i);
    for j = tracker : numel(search_series)
        if (search_series(j)-tmp) <= 0 && (search_series(j+1)-tmp) >0
            if abs (search_series(j)-tmp) < abs (search_series(j+1)-tmp)
                indices(i) = j;
                is_indices(j) = 1;
                tracker = j;
                continue;
            else
                indices(i) = j+1;
                is_indices(j+1) = 1;

                tracker = j+1;
                continue
            end
        end
    end
end
%%
is_indices = is_indices(is_indices == 1);
times_found = search_series(is_indices);
end

