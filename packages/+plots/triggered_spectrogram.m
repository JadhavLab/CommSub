function triggered_spectrogram(efizz, specs, oneDplots, varargin)
% plots triggered spectrogram data from efizz

ip = inputParser;
ip.addParameter('time', [], @isnumeric);
ip.addParameter('thresholds', [], @isstruct);
ip.addParameter('folder', 'commsub_triger', @ischar);
ip.addParameter('zscore', [], @(x) iscellstr(x) || isstring(x));
ip.addParameter('nolog', [], @(x) iscellstr(x) || isstring(x));
ip.parse(varargin{:});
Opt = ip.Results;

folder = figuredefine(Opt.folder);

% figure
fields = fieldnames(specs);
if isempty(Opt.time)
    Opt.time = 1:size(specs.(fields{1}),1);
end
tiledlayout(length(fields)+length(oneDplots), 1)
fg = gcf;
axs = [];
for i = 1:length(fields)
    axs(i) = nexttile;
    f = fields{i};
    if ismember(f, Opt.zscore)
        % along the time axis
        specs.(f) = zscore(specs.(f), [], 1);
        imagesc(Opt.time, efizz.f, specs.(f)')
    elseif ismember(f, Opt.nolog)
        imagesc(Opt.time, efizz.f, specs.(f)')
    else
        imagesc(Opt.time, efizz.f, log(specs.(f))')
    end
    set(gca, 'YDir', 'normal')
    title(f)
end
linkaxes(axs, 'xy')
if ~isempty(oneDplots)
    nexttile
    for i = 1:length(oneDplots)
        if iscell(oneDplots)
            odp = oneDplots{i};
        else
            odp = oneDplots(i);
        end
        fields = fieldnames(odp);
        for j = 1:length(fields)
            if ismember(fields{j}, Opt.zscore)
                % along the time axis
                odp.(fields{j}) = zscore(odp.(fields{j}), [], 1);
            end
            if ~isempty(Opt.thresholds) && isfield(Opt.thresholds, fields{j})
                yline(Opt.thresholds.(fields{j}), 'r')
            end
            plot(Opt.time, odp.(fields{j}))
            hold on
        end
        set(gca,'xlim', [min(Opt.time) max(Opt.time)])
        legend(fields)
        title(strjoin(fields, ' '))
    end
    linkaxes(findobj(fg, 'Type', 'Axes'), 'x')
end
