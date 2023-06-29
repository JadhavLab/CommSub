function [rawLabels] = internal(matrix, labeler_dimension, total_dimension, hotencoder_example, hotencoder_original, ...
    inherent_repelem, round_precision, varargin)
% LABLER.INTERNAL

%% Optional input
ip = inputParser;
ip.addParameter('summary_request', [])
ip.parse(varargin{:});
opt = ip.Results;


%% Determine summary function
if isempty(opt.summary_request)
    if isnumeric(hotencoder_original)
        summary = @mean;
    elseif isstring(hotencoder_original) || ischar(hotencoder_original)
        summary = @mode;
        hotencoder_original = categorical(hotencoder_original);
    elseif iscategorical(hotencoder_original)
        summary = @mode;
    end
else
    summary = opt.summary_request;
end

%% Remove points that will throw off summary
removals = isnan(hotencoder_original) | isinf(hotencoder_original);
removals = removals | logical(sum(isnan(hotencoder_example),2));
hotencoder_original(removals) = [];
hotencoder_example(removals,:) = [];


%% Find the labels from the hotencoder example
if inherent_repelem
    hotencoder_example = hotencoder_example(:, 1:inherent_repelem:end);
else
    disp('no inherent_repelem')
end

assert(all(ismember(nansum(hotencoder_example,2),[0 1])), 'Uh oh, not a pure hotencoding. Please examine hotencoder_example')
[label_hotencoded, ~] = find(hotencoder_example');
assert(any(numel(label_hotencoded) == size(hotencoder_original)), 'Uh oh, dimension mismatch!')
uLabel_hotencoded = unique(label_hotencoded);
% Determine the raw labels from the original sequence represented by the hotencoded matrix?
rosettaStone = zeros(size(uLabel_hotencoded));
for lab = uLabel_hotencoded'
    rosettaStone(lab) = summary(hotencoder_original(lab == label_hotencoded));
end
if round_precision
    rosettaStone = round(rosettaStone, round_precision);
end

%% Finalize the labels
nLabelerDim = size(matrix, labeler_dimension);
assert(nLabelerDim >= numel(uLabel_hotencoded) * inherent_repelem)
rawLabels = rosettaStone;
if inherent_repelem
    rawLabels = repelem(rawLabels, inherent_repelem);
end

