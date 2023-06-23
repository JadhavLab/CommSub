function plotSingleEvents(cellOfWindows, ...
                          Events, Option,  r, varargin)
%   PLOTSINGLEEVENTS(CELLOFWINDOWS, EVENTS, OPTION, R, VARARGIN)
%   Plots the events in CELLOFWINDOWS for the given EVENTTYPE.
% Inputs:
%   cellOfWindows   - Cell array of windows for each event type
%   Events          - Events structure
%   Option          - Option structure
%   r               - Raster structure
%   varargin        - 'displayType' - 'heatmap' or 'raster'
%                   must be followed by the display type
%                   Default == 'heatmap'
%                   'a' - Number of seconds to add to the end of the window
%                   must be a numeric scalar Default == 1
%                   'b' - Number of seconds to subtract from the beginning
%                   of the window must be a numeric scalar Default == 1
%                   'eventType' - Event type to plot
%                   must be a string scalar or a character vector
%                   Default == 'theta'
%   other varargin options are passed to plots.event.processEventsPreamble

ip = inputParser();
ip.KeepUnmatched = true;
ip.addParameter('displayType', 'heatmap', @ischar);
ip.addParameter('after', 1, @isscalar); % Number of seconds to add to the end of the window
ip.addParameter('before', 1, @isscalar); % Number of seconds to subtract from the beginning of the window
ip.addParameter('eventTypes', 'theta')
ip.addParameter('saveFolder', fullfile(figuredefine(), 'singleEvents'),...
                @ischar);
ip.addParameter('appendFigTitle', '', @ischar);
ip.addParameter('subspaceComponents', [], @isnumeric);
ip.parse(varargin{:});
Opt = ip.Results;
varargin = {ip.Unmatched};

if ~exist(Opt.saveFolder, 'dir')
    mkdir(Opt.saveFolder)
end

[Opt2, H, nPatterns, ~, windowQ, windowQlow] = ...
    plots.event.processEventsPreamble(Events, Option, varargin{:});

for field = fieldnames(Opt2)
    Opt.(field{1}) = Opt2.(field{1});
end

for eventType = string(Opt.eventTypes)

    % Determine which event type is considered
    switch eventType
    case "theta"
        eventIdx = 1;
    case "delta"
        eventIdx = 2;
    case "ripple"
        eventIdx = 3;
    otherwise
        error('Unknown event type.')
    end

    % Assume cellOfWindows is a cell array with each cell containing [start, end] 
    % of a window for the chosen eventType
    windows = cellOfWindows{eventIdx};

    disp("Plotting " + num2str(size(windows,1)) + " windows " ...
    + " for " + eventType + " events.");


    % For each window
    for i = 1:size(windows,1)

        % Get the window
        window     = windows(i, :);
        boundary   = [window(1)-Opt.before, window(2)+Opt.after];
        event_inds = find(Events.times > boundary(1) ...
                            & Events.times < boundary(2));

        % Get event times and values
        eventTimes = Events.times(event_inds);
        eventVals  = Events.Hvals(event_inds, eventIdx);

        clf

        % Plot the event values
        subplot(2, 1, 1);
        cla
        plot(eventTimes, eventVals, 'k');
        hold on;
        % Horizontal line for threshold
        plot([min(xlim()), max(xlim())], ...
            [windowQ(eventIdx), windowQ(eventIdx)], 'r');
        % Horizontal line for low threshold
        plot([min(xlim()), max(xlim())], ...
            [windowQlow(eventIdx), windowQlow(eventIdx)], 'r', ...
                 'LineStyle', '--');
        % Plot the window with transparent span
        fill([window(1), window(1), window(2), window(2)],...
             [min(ylim()), max(ylim()), max(ylim()), min(ylim())],...
                'r', 'FaceAlpha', 0.1);
        hold off;

        % Plot the spikes
        ax=subplot(2, 1, 2);
        switch Opt.displayType
            case 'raster'
                error('Raster not implemented yet.')
            case 'heatmap'
                % Here you will plot your spike counts as a heatmap, depends on how you have them stored
                cla(ax)
                spike_count_times  = r.timeBinMidPoints;
                spike_count_inds = find(spike_count_times > boundary(1)...
                            & spike_count_times < boundary(2));
                spikeCountTimes = spike_count_times(spike_count_inds);
                spikeCounts = r.spikeCountMatrix(:, spike_count_inds);  % modify this according to your data
                plots.heatmapBinSpikes(spikeCounts, 'times', spikeCountTimes,...
                    'cmap', flipud(cmocean('solar')), 'colorbar', true, ...
                                       'background', 'white', 'ax', ax);
                hold on;
                fi=fill(ax,[window(1), window(1), window(2), window(2)],...
                    [min(ylim()), max(ylim()), max(ylim()), min(ylim())], ...
                    'r', 'FaceAlpha', 0.1);
                uistack(fi,'top');

            otherwise
                hold on;
                error('Unknown display type.')
        end % switch Opt.displayType

        if ~isempty(Opt.subspaceComponents)
            keyboard;
            spike_count_times  = r.timeBinMidPoints;
            spike_count_inds = find(spike_count_times > boundary(1)...
                        & spike_count_times < boundary(2));
            spikeCountTimes = spike_count_times(spike_count_inds);
            objs = []
            for j = 1:size(Opt.subspaceComponents, 1)
                sub = Opt.subspaceComponents(j, event_inds)
                sub = (sub - min(sub)) / (max(sub) - min(sub)) .* ...
                    (max(ylim()) - min(ylim())) + min(ylim())
                subplot(2, 1, 1);
                hold on;
                o1=plot(spike_count_times, sub)
                    
                subplot(2, 1, 2);
                hold on;
                o2=plot(spike_count_times, sub)
                objs = [objs, o1, o2];
            end
        end

        linkaxes([subplot(2, 1, 1), subplot(2, 1, 2)], 'x');
        % aspect ratio 0.3
        set(gcf, 'Position', [100, 100, 300, 900]);
        sgtitle(eventType + " event " + num2str(i) + " of " + ...
            num2str(size(windows,1)) + newline + Opt.appendFigTitle);

        % Save the figure
        saveas(gcf, fullfile(Opt.saveFolder, ...
            eventType + "_" + num2str(i) + ".png"));
        saveas(gcf, fullfile(Opt.saveFolder, ...
            eventType + "_" + num2str(i) + ".svg"));

    end % for each window

end % for each event type
