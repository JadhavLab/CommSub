function eventTypes = addEventType(cellOfWindows, animal_behavior)

% this function takes in the animal's behavior/component time-wise table
% and output the type of event that is happening during each of the window

% v1: linear search from the beginning every time, just a brute-force
% solution
eventTypes = zeros(1,height(animal_behavior));
for i = 1:height(animal_behavior)
    currTime = animal_behavior.time(i);
    for j = 1:numel(cellOfWindows)
        currSetOfWindows = cellOfWindows{j};
        for k = 1:size(currSetOfWindows,1)
            if currTime >= currSetOfWindows(k,1) & currTime < currSetOfWindows(k,2)
                eventTypes(i) = j;
                continue
            end
        end
    end
end

end

