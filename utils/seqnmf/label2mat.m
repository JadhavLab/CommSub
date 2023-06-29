function mat = label2mat(label, size, varargin)
% label2mat one hot encoder
% 
% label : time series to hot encode
% size  : size to hot encode the time series into 

ip = inputParser;
ip.addParameter('method', []);
ip.parse(varargin{:});
opt = ip.Results;

% Determine size if not given?
if ~exist('size', 'var') || isempty(size)
    size = numel(unique(label));
end

% Discretize?
if isempty(opt.method)
    isReal = ~nansum( round(label(:)) - label(:) ) == 0;
    if isReal
        label = discretize(label,size);
    end
else
    switch opt.method
    case 'discretize'
        label = discretize(label,size);
    case 'quantile'
        quantile_edges = quantile(label, (1:(size+1))/(size+1));
        label = discretize(label, quantile_edges);
    case 'raw'
        % Label is already dicretize into bins
    otherwise
        error("Unrecognized method= " + opt.method)
    end
end

if label > size
    error('Label (%d) should be < size (%d).', label, size);
end

% Create hot encoding
I = eye(size);
nanlocs = isnan(label);
label(nanlocs) = 1;
mat = I(:, label);
mat(:,nanlocs) = nan;

