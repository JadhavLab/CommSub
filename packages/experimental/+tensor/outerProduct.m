function [z, times] = outerProduct(x, y, times, varargin)

ip = inputParser;
ip.addParameter('center',true);
ip.addParameter('downsample', 0);
ip.addParameter('cast', "single");
ip.addParameter('timedim',-1);
ip.addParameter('nontimedim',-1);
p.parse(varargin{:})
opt = ip.Results;

if opt.timedim == -1
    opt.timedim(1) = ndims(x);
    opt.timedim(2) = ndims(y);
end
if opt.nontimedim == -1 || ~iscell(opt.nontimedim)
    opt.nontimedim = {};
    opt.nontimedim{1} = setdiff(1:ndims(x), opt.timedim(1));
    opt.nontimedim{2} = setdiff(1:ndims(x), opt.timedim(2));
    opt.nontimedim = cellfun(@(x) num2cell(x), opt.nontimedim, 'UniformOutput', false);
end
nontimesize{1} = size(x, [opt.nontimedim{1}{:}]);
nontimesize{2} = size(y, [opt.nontimedim{2}{:}]);

if opt.center
    disp('centering')
    x = center(x, opt.timedim(1));
    y = center(y, opt.timedim(2));
end

disp('tensor to matrix')
x = tens2mat(x, opt.timedim(1));
y = tens2mat(y, opt.timedim(2));
if opt.cast == "single"
    x = single(x);
    y = single(y);
end
if opt.downsample ~= 0
    %x = downsample(x, opt.downsample);
    %y = downsample(y, opt.downsample);
    %times = downsample(times, opt.downsample);
    %opt.dow5nsample = timelen(x,1)/opt.downsample;
    if opt.downsample > 0
        kws = {'Subsample', opt.downsample};
    else
        kws = {};
    end
    x = decimate(x, kws{:}, 'Full', true);
    x = squeeze(mean(x,2));
    %x = mean(x, ndims(x));
    y = decimate(y, kws{:}, 'Full', true);
    y = squeeze(mean(y,2));
    %y = mean(y, ndims(y));
    times = decimate(times, kws{:}, 'Full', true);
    times = squeeze(mean(times,2));
    %times = mean(times, ndims(times));
end
timelen = size(x,1);
y = shiftdim(y, -1);
y = permute(y, [2 1 3]);
z = bsxfun(@times, x, y);
z = reshape(z, [timelen, nontimesize{1}, nontimesize{2}]);



function x = center(x, timedim);

    x = (x - nanmean(x, timedim))./nanstd(x, 1, timedim);
    x(isinf(x)) = 0;
