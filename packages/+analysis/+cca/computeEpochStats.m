function statsTable = computeEpochStats(cca, behavior, Spk)
    % Initialize a table to store the results
    statsTable = table();
    
    % Interpolate U and V to match behavior timestamps
    cca.U = interp1(cca.times, cca.U, behavior.times);
    cca.V = interp1(cca.times, cca.V, behavior.times);
    
    % Interpolate spikeRateMatrix to match behavior timestamps
    spikeRateMatrix = interp1(Spk.timeBinMidPoints, Spk.spikeRateMatrix, behavior.times);
    
    % Find the indices where trajbound changes
    changePoints = find(diff(behavior.trajbound) ~= 0);
    
    % Add the start and end points
    changePoints = [0; changePoints; length(behavior.trajbound)];
    
    % For each epoch, compute the desired statistics
    for i = 1:(length(changePoints) - 1)
        % Extract the epoch
        startIdx = changePoints(i) + 1;
        endIdx = changePoints(i + 1);
        epochBehavior = behavior(startIdx:endIdx, :);
        epochU = cca.U(:, startIdx:endIdx);
        epochV = cca.V(:, startIdx:endIdx);
        epochMUA = mean(spikeRateMatrix(:, startIdx:endIdx), 1);
        
        % Normalize U and V by MUA
        epochU = epochU ./ epochMUA;
        epochV = epochV ./ epochMUA;
        
        % Compute the statistics
        meanU = mean(epochU, 2);
        medianU = median(epochU, 2);
        meanV = mean(epochV, 2);
        medianV = median(epochV, 2);
        
        % Create a row in the results table
        resultRow = table(epochBehavior.trajbound(1), meanU, medianU, meanV, medianV, ...
                          'VariableNames', {'Epoch', 'MeanU', 'MedianU', 'MeanV', 'MedianV'});
        
        % Append the result row to the results table
        statsTable = [statsTable; resultRow];
    end
end

