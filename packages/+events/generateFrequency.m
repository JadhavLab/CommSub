function frequency = generateFrequency(animal, type_of_frequency, brainArea, varargin)


% this function generate the event matrix for the specified osciallation
% matrix


if type_of_frequency == "ripple"
    load(animal + "globalripple01.mat");
    
      [~, H(:,RIPPLE), Hnanlocs(:,RIPPLE), Hvals(:,RIPPLE), minRippleThreshold, original] = ...
            eventMatrix.generateFromRipples(globalripple, ...
            'amplitude_at_riptime', true,...
            'rippleBand', Hvals(:,RIPPLE),... RY: Hvals, not H here for obvious reasons: you  want the original ripple band activity
            'rippleBandTime', Htimes,...
            'globalrippleWindowUnits', 'std');
    % do the times from 
    
    
else
    load(animal + "avgeeg.mat");
    % only take the running session data
    [runningSess, ~] = getRunningSessions(animal);
    % If user passes a branched cell array, instead of an ndimension struct, converet it!
    if iscell(avgEEGStruct)
        indices = ndBranch.indicesMatrixForm(avgEEGStruct); %get the indices of all the cells within the avg eeg files, each with pfc and hpc info
        avgEEGStruct = ndBranch.toNd(avgEEGStruct); %convert into structs with relevant information
    else
        indices = nd.indicesMatrixForm(avgEEGStruct);
    end
    
    for index = indices'
        if ~ismember(index(2), runningSess)
            continue
        end
        
        I = num2cell(index);
        singleStruct = avgEEGStruct( I{:} );
        
        % If this data not from the requested brain area, skip
        if ~isequal( singleStruct.area , brainArea )
            continue
        end
        
        % Get times of this single day-ep-brainarea
        singleTime = geteegtimes(singleStruct);
        % Data for ths single  day-ep-brainarea
        eegData     = singleStruct.(type_of_frequency).data;
        eegData = double(eegData);
        
        frequency = [singleTime',eegData(:,2)/10000,eegData(:,1)];       
    end
end

end

