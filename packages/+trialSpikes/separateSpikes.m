function X_regionPattern = separateSpikes(spikeSample, areaPerNeuron, region)
% this function separate the spikeSample of each neuron based on its specified
% region, along with the corresponding firing rate of these neurons
%
% spikeSample can be a matrix or a tensor

nPatterns = numel(spikeSample);
X_regionPattern = cell(size(spikeSample));

%% for each pattern, get the neurons
for r = 1:nPatterns
    curr = spikeSample{r};
    if ismatrix(curr)
        % [nNeurons,~] = size(curr);
    elseif ndims(curr) == 3
        % [nNeurons,~, ~] = size(curr);
    else
        error("Unsupported dimension of spikeSample");
    end
    currResult = curr(areaPerNeuron == region, :, :);
    X_regionPattern{r} = currResult;
end
