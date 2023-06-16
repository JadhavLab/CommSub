function [equalized, warned] = equalizePatternControl(cellOfWindows)
% equalized the number of windows between a pattern and its control
% cellOfWindows has to be formatted so that the first nPatterns are the
% patterns and the second half are the controls

equalized = cellOfWindows;
warned = false;

for pattern = 1:length(cellOfWindows)/2
    curr_pattern = cell2mat(cellOfWindows(:,pattern));
    curr_control = cell2mat(cellOfWindows(:,pattern+length(cellOfWindows)/2));
    
    % if control window pattern is empty (due to removing overlaps), use
    % the original number of the windows in the pattern.
    if isempty(curr_control)
        warning("control" + pattern +...
                " is empty,using original number of pattern windows");
        warned = true;
        break;
    end
    
    [real_pattern, real_control] = ...
        windows.removeUntilEqual(curr_pattern, curr_control, "matched", 50);
    equalized{pattern} = real_pattern;
    equalized{pattern+length(cellOfWindows)/2} = real_control;
end


end

