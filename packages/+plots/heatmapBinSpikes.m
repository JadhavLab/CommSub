function heatmapBinSpikes(r, varargin)
% heatmapBinSpikes(r, varargin)
%
% Plots a heatmap of the spike count matrix.
%
% Inputs:
%   r           - the results struct | spikeCountMatrix
%
% Options:
%   ax          - axis to plot on
%   colorCols   - which color columns to use for the heatmap

    ip = inputParser();
    ip.addParameter('ax', gca, @ishandle);
    ip.addParameter('colorCols', 1:3, @isnumeric);
    ip.addParameter('ylim', [], @isnumeric);
    ip.addParameter('times', [], @isnumeric);
    ip.addParameter('cmap', gray(256), @isnumeric);
    ip.addParameter('colorbar', true, @islogical);
    ip.addParameter('upperHalfCmap', false, @islogical);
    ip.addParameter('inverseCmap', false, @islogical);
    ip.addParameter('background', 'w', @(x) ischar(x) || isnumeric(x));
    ip.parse(varargin{:});
    Opt = ip.Results;

    if isstruct(r)
        scm   = r.spikeCountMatrix;
    else
        scm   = r;
    end
    if isempty(Opt.times) && isstruct(r)
        times = r.timeBinMidPoints;
    elseif isempty(Opt.times)
        times = 1:size(scm, 2);
    else
        times = Opt.times;
    end

    % Compute the color map
    cmap = Opt.cmap;
    % cmap = flipud(cmap);
    cmap(:, setdiff(1:3,Opt.colorCols)) = 0;
    if Opt.upperHalfCmap
        cmap = cmap(128:end, :);
        cmap = repelem(cmap, 2, 1);
    end
    if Opt.inverseCmap
        cmap = flipud(cmap);
    end
    cmap(1, :) = 1;

    % ylimits
    ylim_given = false;
    if ~isempty(Opt.ylim) 
        ylim_given = true;
        if numel(Opt.ylim) == 2
            Opt.ylim = linspace(Opt.ylim(1), Opt.ylim(2), size(scm, 1));
        end
    else
        Opt.ylim = 1:size(scm, 1);
    end

    % Plot the heatmap
    im = imagesc(Opt.ax, times, Opt.ylim, scm);
    if ~ylim_given
        ylim([min(Opt.ylim), max(Opt.ylim)]);
        ylabel("Cell #");
    end

    uistack(im, 'bottom');
    colormap(Opt.ax, cmap);
    % caxis([0 1]);
    if ischar(Opt.background)
        [name,Opt.background] = colornames("wikipedia", Opt.background);
        disp("SEtting background to " + name);
    end
    set(gca, 'YDir', 'normal', 'Color', Opt.background);
    % xlim([times(1) times(end)]);
    if Opt.colorbar
        colorbar('off');
        cb = colorbar('Location', 'eastoutside');
        cb.Label.String = 'Spike Count';
    end
