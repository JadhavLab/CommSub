function triggered_spectrogram(efizz, specs, oneDplots, varargin)
% plots triggered spectrogram data from efizz

ip = inputParser;
ip.addParameter('specNames', {'S1', 'S2', 'C'}, @iscell);
ip.addParameter('time', 1:size(specs{1},1), @isnumeric);
ip.addParameter('oneDleg', {{'u', 'v'}}, @iscell);
ip.addParameter('oneDnames', {{'u and v'}}, @iscell);
ip.addParajjkjkjkmeter('folder', 'commsub_triger', @ischar);
ip.add
ip.parse(varargin{:});
Opt = ip.Results;

folder = figuredefine(Opt.folder);

tiledlayout(length(Opt.specNames)+~isempty(oneDplots), 1)
for i = 1:length(Opt.specNames)
    nexttile
    imagesc(1:size(specs{i},1), efizz.f, log(specs{i}))
    set(gca, 'YDir', 'normal')
    title(Opt.specNames{i})
end
if ~isempty(oneDplots)
    nexttile
    for i = 1:length(oneDplots)
        plot(Opt.time, oneDplots{i})
        hold on
    end
    legend(Opt.oneDleg{i}{:})
    title(Opt.oneDnames{i})
end
