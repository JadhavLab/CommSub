function [newTimes, animalStartTimes]=makeSuperAnimal(times)
    %OVERALL FUNCTION PURPOSE
    % 
    %OUTPUTS
    %newTimes: a vector of times points with several animal's
    %          experimental time merged
    %animalStartTimes: start time for each animal : list of times in newTimes 
    % that are the first time of new animal data 
    % (helps you know which times actually are new animals).
    %INPUT
    %times: time vector 
    
    offset = 100000; % Amount of time to offset each new animal by
    
    D = diff(times);

    animalStartTimes = find(D<0);
    D(D<0) = offset;
    newTimes = cumsum(D);
end

% this is needed because for every recording the time starts at 0