function [newSpikeCountMatrix, newTimeBinMidPoints] = fillSleepPeriods(spikeCountMatrix, timeBinMidPoints)


    ip = inputParser();
    ip.addParameter('medianSampleRateFactor', 20);
    ip.parse(varargin{:});
    Opt = ip.Results;


    % Calculate the sample rate
    timeDiffs = diff(timeBinMidPoints);
    medianSampleRate = median(timeDiffs);

    % Identify the sleep periods
    sleepPeriods = find(timeDiffs > (medianSampleRate * Opt.medianSampleRateFactor));
    
    % Start with the original data
    newSpikeCountMatrix = spikeCountMatrix;
    newTimeBinMidPoints = timeBinMidPoints;

    % Loop backwards through the sleep periods
    for i = flip(sleepPeriods)
        % Calculate the number of missing samples
        missingSamples = round(timeDiffs(i) / medianSampleRate) - 1;

        % Add NaN values for the missing samples in the spike count matrix
        newSpikeCountMatrix = [newSpikeCountMatrix(:, 1:i), nan(size(spikeCountMatrix, 1), missingSamples), newSpikeCountMatrix(:, i+1:end)];

        % Add the missing times in the time bin midpoints
        missingTimes = timeBinMidPoints(i) + medianSampleRate * (1:missingSamples)';
        newTimeBinMidPoints = [newTimeBinMidPoints(1:i); missingTimes; newTimeBinMidPoints(i+1:end)];
    end

end

