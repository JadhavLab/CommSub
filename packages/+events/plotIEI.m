function plotIEI(Events, varargin)
    ip = inputParser;
    ip.addParameter('quantile', 0.01, @isnumeric); % quantile to throw out
    ip.addParameter('layout', 'single', @(x) ismember(x, {'single', 'row'})); % layout of histograms
    ip.KeepUnmatched = true;
    ip.parse(varargin{:});
    Opt = ip.Results;
    varargin = [fieldnames(ip.Unmatched) struct2cell(ip.Unmatched)]';

    % Get the number of event sets
    numEventSets = length(Events.cellOfWindows);

    if strcmp(Opt.layout, 'single')
        % Create a new figure
        figure;
        hold on;

        % Iterate over each event set
        for i = 1:numEventSets
            % Generate and plot the histogram for the current event set
            plotIEIHistogram(Events, i, Opt, varargin);
        end

        hold off;
    elseif strcmp(Opt.layout, 'row')
        % Create a new figure with a row of subplots
        figure;
        for i = 1:numEventSets
            % Create a subplot for the current event set
            subplot(1, numEventSets, i);
            hold on;

            % Generate and plot the histogram for the current event set
            plotIEIHistogram(Events, i, Opt, varargin);

            hold off;
        end
    end

    % Add labels and a legend to the figure
    xlabel('IEI');
    ylabel('Probability');
    legend(Events.cellOfWin_varnames);
end

function plotIEIHistogram(Events, i, Opt, varargin)
    % Get the current set of events
    currEvents = Events.cellOfWindows{i};

    % Compute the start times of the events
    startTimes = currEvents(:, 1);

    % Compute the IEI
    IEI = diff(startTimes);

    % Throw out the top and bottom quantile
    IEI = IEI(IEI > quantile(IEI, Opt.quantile) & IEI < quantile(IEI, 1 - Opt.quantile));

    % Generate a histogram of the IEI
    histogram(IEI, 'Normalization', 'probability', varargin{:});

    % Add a title to the histogram
    title(sprintf('IEI Histogram for %s', Events.cellOfWin_varnames{i}));
end

