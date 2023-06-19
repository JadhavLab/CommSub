function heatmapBinSpikes(r, varargin)
% heatmapBinSpikes(r, varargin)
%
% Plots a heatmap of the spike count matrix.
%
% Inputs:
%   r           - the results struct
%
% Options:
%   ax          - axis to plot on
%   colorCols   - which color columns to use for the heatmap

    ip = inputParser();
    ip.addParameter('ax', gca, @ishandle);
    ip.addParameter('colorCols', 1:3, @isnumeric);
    ip.addParameter('ylim', [], @isnumeric);
    ip.parse(varargin{:});
    Opt = ip.Results;

    scm   = r.spikeCountMatrix;
    times = r.timeBinMidPoints;

    % Compute the color map
    cmap = gray;
    % cmap = flipud(cmap);
    cmap(:, setdiff(1:3,Opt.colorCols)) = 0;
    cmap = cmap(128:end, :);
    cmap = repelem(cmap, 2, 1);
    cmap(1, :) = 1;

    % Plot the heatmap
    axes(Opt.ax);
    % hold on;
    if isempty(Opt.ylim)
        im=imagesc( Opt.ax, times, 1:height(r.celllookup), r.spikeCountMatrix);
        ylim([1 height(r.celllookup)]);
        ylabel("Cell #");
    else
        im=imagesc( Opt.ax, times, Opt.ylim, r.spikeCountMatrix);
    end
    uistack(im, 'bottom');
    colormap(Opt.ax, cmap);
    % caxis([0 1]);
    set(gca, 'YDir', 'normal');
    xlim([times(1) times(end)]);
