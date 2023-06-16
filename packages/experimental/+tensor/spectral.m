function [values, times] = lfp(efizz, hankel, varargin)
% Returns a hankelized tensor of lfp

ip = inputParser;
ip.addParameter('fields',["S1","S2","C"]);
ip.addParameter('centerTimes', false);
ip.addParameter('query_time', []);
ip.addParameter('squeeze',true);
ip.addParameter('single',true);
ip.addParameter('center',true);
ip.parse(varargin{:})
Opt = ip.Results;

times = efizz.t;
values = cell(numel(Opt.fields),1);
for field = string(Opt.fields)

    if startsWith(field,'S')
        efizz.(field) = abs(efizz.(field));
    end
    if Opt.center
        values{end+1} = efizz.(field);
    else
        values{end+1} = center(efizz.(field), 1);
    end
end
values = cat(2, values{:});
values = values';

if  ~isempty(Opt.query_time)
    inds = interp1(times, 1:numel(times), Opt.query_time);
    values = values(inds,:);
    times = times(inds);
end

if ~isempty(hankel)
    values = hankelize(values, 'Sizes', hankel, 'Full', true, 'Dim', 2);
    times =  hankelize(times, 'Sizes',   hankel, 'Full', true);
end

if Opt.squeeze
    if iscell(values)
        values = cellfun(@(value) squeeze(value), values);
    else
        values = squeeze(values);
    end
end

if Opt.single
    if iscell(values)
        values = cellfun(@(value) single(value), values);
    else
        values = single(values);
    end
end


if Opt.centerTimes
    times = median(times,2);
end

function x = center(x, timedim);

    x = (x - nanmean(x, timedim))./nanstd(x, 1, timedim);
    x(isinf(x)) = 0;
