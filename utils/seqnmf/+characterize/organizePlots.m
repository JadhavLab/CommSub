function organizePlots(data, func, varargin)
% ORGANIZEPLOTS allows user to organize either a wholistic, K separated or F separated plot

ip = inputParser;
ip.addParameter('vartype', 'W'); % W | WH | '' | H
ip.addParameter('wl_kws', {'downsample', 30, 'protectR', [0, 1, -1], 'protectField', {'velocity'}); % Default KWS for wholeLabel
ip.addParameter('collapseF', true); % collapse all P
ip.addParameter('collapseK', true); % collapse all K 
ip.addParameter('colorscheme', 'hspace'); % Default KWS for wholeLabel
ip.addParameter('colorscheme', 'vspace'); % Default KWS for wholeLabel
ip.parse(varargin{:});
opt = ip.Results;



% Get the dimension that corresponds to time and K
td = timedim(vartype, data.fields{1});
kd = kdim(vartype);

% Obtain labels for axis
[axP, axR, axisCenters] = labeler.wholeLabel(data, data.fields, opt.wholeLabel_kws{:}, 'collapse', opt.collapseF);

%% Axis shortcuts
% Determine the axis to spread across for the k dimension
if opt.kAxis == 1
    nestK = @(K, k) nestplot(K, 1, k);
    nest = @(K, F, k, f) nestplot(K, F, k, f);
elseif opt.kAxis == 2
    nestK = @(K, x) nestplot(1, K, x);
    nest = @(K, F, k, f) nestplot(F, K, f, k);
end
% Determine the axis to spread across for the k dimension
if opt.fAxis == 1
    nestK = @(F, f) nestplot(F, 1, f);
elseif opt.fAxis == 2
    nestK = @(F, f) nestplot(1, F, f);
end

%% Setup playfield
axFKtracker = {};
if opt.collapseF && opt.collapseK

    axNest = gca;
    axNest = nestable(axNest);
    F = numel(data.fields);
    K = size(data.W,2);
    for f = 1:F
        for k = 1:K
            axFFtracker{k,f} = nest(K, F, k, f);
            % Pull out the field
            D = data.(data.field{F})
        end
    end

elseif opt.collapseF

    axNest = gca;
    axNest = nestable(axNest);
    F = numel(data.fields);
    for k = 1:K
        for f = 1:F
            axFFtracker{1,f} = nest(K, F, k, f);
        end
    end

elseif opt.collapseK

    axNest = gca;
    axNest = nestable(axNest);
    K = size(data.W,2);
    for k = 1:K
        axFKtracker{k,1} = nest(K, k);
    end

else
    varargin = [varargin, 'axisCenters',axisCenters, 'axR', axR, 'axP', axP];
    func(data, varargin);
end


