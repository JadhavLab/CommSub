function [first_time_of_the_day, last_time_of_the_day] = getTimerange(ripples)

% this function gets the first and last timepoint of a eeg file
% ripples is an array of ripple records in the format of (start, stop, amp)
first_time_of_the_day = intmax;
last_time_of_the_day = -1;

% find the inital and end time of the day
for iEpoch = 1:numel(ripples)
    if ~isempty(ripples{iEpoch})
        curr = ripples{iEpoch};
        curr_first = curr(1,1);
        curr_last = curr(end,2);
        if curr_first<first_time_of_the_day
            first_time_of_the_day = curr_first;
        end
        if curr_last>last_time_of_the_day
            last_time_of_the_day = curr_last;
        end
    end
end

end

