function plotEpochStats(statsTable)
    % Get the number of unique epochs
    numEpochs = numel(unique(statsTable.Epoch));
    
    % Get the number of U and V components
    numComponents = size(statsTable.MeanU, 1);
    
    % Create a new figure
    figure;
    
    % For each component, create a subplot and plot the statistics
    for i = 1:numComponents
        % Create a subplot for U
        subplot(numComponents, 2, 2*i-1);
        % Plot mean and median U
        bar([statsTable.MeanU(i, :); statsTable.MedianU(i, :)]);
        % Set title and labels
        title(sprintf('U component %d', i));
        xlabel('Epoch');
        ylabel('Value');
        legend('Mean', 'Median');
        
        % Create a subplot for V
        subplot(numComponents, 2, 2*i);
        % Plot mean and median V
        bar([statsTable.MeanV(i, :); statsTable.MedianV(i, :)]);
        % Set title and labels
        title(sprintf('V component %d', i));
        xlabel('Epoch');
        ylabel('Value');
        legend('Mean', 'Median');
    end
end

