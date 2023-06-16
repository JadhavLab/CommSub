function [values, times] = spikes(values, times, hankel, varargin)
% Returns a hankelized tensor of spikes

ip = inputParser;
ip.addParameter('areaPerNeuron',[]);
ip.addParameter('query_time',[]);
ip.addParameter('Dim',1); %Dimension to hankelize over
ip.addParameter('centerTimes',false);
ip.addParameter('squeeze',true);
ip.addParameter('single',true);
ip.parse(varargin{:})
Opt = ip.Results;

% Get neurons by area? one cell per area
if ~isempty(Opt.areaPerNeuron)
    V = values;
    uVals = unique(Opt.areaPerNeuron);
    values = cell(numel(uVals),1);
    for uVal = uVals
        values{i} = V(Opt.areaPerNeuron == uVal,:);
    end
end

% Interpolate to some passed times?
if nargin > 3 && ~isempty(Opt.query_time)
    inds = interp1(times, 1:numel(times), Opt.query_time);
    values = values(inds,:);
    times = times(inds);
end

% Hanelize
if ~isempty(hankel)
    if iscell(values)
        values = cellfun(@(value) hankelize(value, 'Sizes', hankel, 'Dim', Opt.Dim), values);
    else
        values = hankelize(values, 'Sizes', hankel, 'Full', true, 'Dim', Opt.Dim);
    end
    times = hankelize(times, 'Sizes', hankel, 'Full', true);
end

%Squeeze?
if Opt.squeeze
    if iscell(values)
        values = cellfun(@(value) squeeze(value), values);
    else
        values = squeeze(values);
    end
end

% Convert to single?
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
