function plotMuaSC(r, varargin)
% plotFR(r, varargin)
% plot firing rate of a neuron

% Plots a heatmap of the spike count matrix.
%
% Inputs:
%   r           - the results struct
%
% Options:
%   ax          - axis to plot on
%   colorCols   - which color columns to use for the heatmap

ip = inputParser();
ip.addParameter('ax',        gca, @ishandle);
ip.addParameter('colorCols', 1:3, @isnumeric);
ip.addParameter('ylim',      [],  @isnumeric);
ip.addParameter('kws',       {},  @iscell);
ip.parse(varargin{:});
Opt = ip.Results;

scm   = sum(r.spikeCountMatrix, 1);
times = r.timeBinMidPoints;

% Plot the heatmap
axes(Opt.ax);
% hold on;
if isempty(Opt.ylim)
    im=plot(times, scm, Opt.kws{:});
else
    scm = (scm - min(Opt.ylim)/(max(Opt.ylim)-min(Opt.ylim))) .* ...
        (max(Opt.ylim)-min(Opt.ylim));
    im=plot(times, scm, Opt.kws{:});
    ylim(Opt.ylim);
end
uistack(im, 'top');
% colormap(Opt.ax, cmap);
% caxis([0 1]);
% set(gca, 'YDir', 'normal');
% xlim([times(1) times(end)]);
