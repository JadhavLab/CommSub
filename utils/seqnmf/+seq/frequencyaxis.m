function frequencyaxis(data, n)

if nargin == 1
    n = 5;
end

freq = data.f(1,:);
spacing = median(diff(freq));
yts = linspace(freq(1)-spacing/2, freq(end)+spacing/2, n)

ytlabels = string(yts);
ytlabels = arrayfun(@(x) sprintf('%2.1f', x), ytlabels, 'UniformOutput', false);
ytlabels = string(ytlabels)

yticks(yts)
yticklabels(ytlabels)
