function plotSingleEvents(eventType, b, a, displayType, cellOfWindows, Events, r)
    % eventType: one of 'theta', 'delta', 'ripple'
    % b: seconds before event
    % a: seconds after event
    % displayType: 'raster' or 'heatmap'
    % cellOfWindows: cell array containing information about time windows around events
    % Events: struct containing information about events
    % r: struct containing information about activity

    eventType = 'theta';
    b = 1;
    a = 1;
    displayType = 'raster';
    
    % Determine which event type is considered
    switch eventType
        case 'theta'
            eventIdx = 1;
        case 'delta'
            eventIdx = 2;
        case 'ripple'
            eventIdx = 3;
        otherwise
            error('Unknown event type.')
    end

    % Assume cellOfWindows is a cell array with each cell containing [start, end] of a window for the chosen eventType
    windows = cellOfWindows{eventIdx};

    % For each window
    for i = 1:size(windows,1)

        % Get the window
        window = windows(i, :);
        boundary = [window(1)-b, window(2)+a];
        event_inds = find(Events.times > boundary(1) & Events.times < boundary(2));

        % Get event times and values
        eventTimes = Events.times(event_inds);
        eventVals  = Events.Hvals(event_inds, eventIdx);

        % Plot the event values
        subplot(2, 1, 1);
        plot(eventTimes, eventVals, 'k');
        hold on;
        % Plot the window with transparent span
        fill([window(1), window(1), window(2), window(2)],...
             [min(ylim()), max(ylim()), max(ylim()), min(ylim())],...
                'r', 'FaceAlpha', 0.1);
        hold off;

        % Plot the spikes
        ax=subplot(2, 1, 2);
        switch displayType
            case 'raster'
                % Here you will plot your spike times, depends on how you have them stored
                % Assume spikeTimes is a matrix where each row corresponds to a neuron and columns are spike times
                spikeTimes = r.spikeRateMatrix;  % modify this according to your data
                for j = 1:size(spikeTimes, 1)
                    plot(spikeTimes(j, :), ones(size(spikeTimes(j, :))) * j, '.k');
                    hold on;
                end
                hold off;

            case 'heatmap'
                % Here you will plot your spike counts as a heatmap, depends on how you have them stored
                cla()
                spike_count_times  = r.timeBinMidPoints;
                spike_count_inds = find(spike_count_times > boundary(1)...
                            & spike_count_times < boundary(2));
                spikeCountTimes = spike_count_times(spike_count_inds);
                spikeCounts = r.spikeCountMatrix(:, spike_count_inds);  % modify this according to your data
                plots.heatmapBinSpikes(spikeCounts, 'times', spikeCountTimes,...
                    'cmap', (cmocean('speed')), 'colorbar', true, 'background', 'black');
                hold on;
                fi=fill(ax,[window(1), window(1), window(2), window(2)],...
                    [min(ylim()), max(ylim()), max(ylim()), min(ylim())], ...
                    'r', 'FaceAlpha', 0.1);
                uistack(fi,'top');

            otherwise
                hold on;
                error('Unknown display type.')
        end

        linkaxes([subplot(2, 1, 1), subplot(2, 1, 2)], 'x');
    end
end

