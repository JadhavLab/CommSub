function plotCrossIEI(Events, varargin)
    ip = inputParser;
    ip.addParameter('quantile', 0.01, @isnumeric); % quantile to throw out
    ip.addParameter('layout', 'matrix', @(x) ismember(x, {'matrix', 'row'})); % layout of histograms
    ip.KeepUnmatched = true;
    ip.parse(varargin{:});
    Opt = ip.Results;
    varargin = [fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]';
    
    % Get the number of event sets
    numEventSets = length(Events.cellOfWindows);
    
    if strcmp(Opt.layout, 'matrix')
        % Create a new figure with a matrix of subplots
        figure;
        for i = 1:numEventSets
            for j = i:numEventSets
                % Create a subplot for the current pair of event sets
                subplot(numEventSets, numEventSets, (i - 1) * numEventSets + j);
                hold on;

                % Generate and plot the histogram for the current pair of event sets
                plotIEIHistogram(Events, i, j, Opt);

                hold off;
            end
        end
    elseif strcmp(Opt.layout, 'row')
        % Create a new figure with a row of subplots
        figure;
        for i = 1:numEventSets
            for j = i:numEventSets
                % Create a subplot for the current pair of event sets
                subplot(1, numEventSets, (i - 1) * numEventSets + j);
                hold on;

                % Generate and plot the histogram for the current pair of event sets
                plotIEIHistogram(Events, i, j, Opt);

                hold off;
            end
        end
    end

    % Add labels and a legend to the figure
    xlabel('IEI');
    ylabel('Probability');
    legend(Events.cellOfWin_varnames);
    linkaxes(findobj(gcf, 'Type', 'Axes'), 'xy');
end

function plotIEIHistogram(Events, i, j, Opt)
    % Get the current pair of event sets
    currEvents1 = Events.cellOfWindows{i};
    currEvents2 = Events.cellOfWindows{j};

    % Compute the start times of the events
    startTimes1 = currEvents1(:, 1);
    startTimes2 = currEvents2(:, 1);

    % Compute the IEI
    IEI = [];
    for k = 1:length(startTimes1)
        % Find the time until the next event in the second set
        timeUntilNextEvent = min(startTimes2(startTimes2 > startTimes1(k))) - startTimes1(k);

        if ~isempty(timeUntilNextEvent)
            IEI = [IEI; timeUntilNextEvent];  % Append to the IEI array
        end
    end

    % Throw out the top and bottom quantiles
    IEI = IEI(IEI > quantile(IEI, Opt.quantile) & IEI < quantile(IEI, 1 - Opt.quantile));

    % Generate a histogram of the IEI
    histogram(IEI, 'Normalization', 'probability', 'NumBins', 100);

    % Add a title to the histogram
    title(sprintf('IEI Histogram for %s and %s', Events.cellOfWin_varnames{i}, Events.cellOfWin_varnames{j}));
end

