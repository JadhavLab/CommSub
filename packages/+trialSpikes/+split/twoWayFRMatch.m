function [Y_source, Y_target, nSource, nTarget, index_source, index_target] =...
    matchFRinDiscreteRanges(X_source, X_target, source_FR, target_FR, nBins, varargin)
%
% this function matches the firing rates of two target region with discrete
% bins from the FR histogram for every pattern
% 
% it will make the prediction matrix more populated by having two regions
% with neurons that have matching firing rate
%
% ------------------------------------------------------------------------%
% input: the firing pattern matrices from the two regions                 %
% output: the matched firing pattern matrices from the two target regions %
%         the number of source and target population eventually decided   %
%         and the neuronal indices of those used                          %
% ------------------------------------------------------------------------%

ip = inputParser;
ip.addParameter('verbose', false, @islogical);
ip.addParameter('ploton',  true, @islogical);
ip.parse(varargin{:});
Opt = ip.Results;


%% (1) sort the firing pattern matrices by firing rate so as to index neurons
% source_FR = hpcFR; target_FR = pfcFR;
numsource = size(source_FR, 2); % number of neurons in source region
numtarget = size(target_FR, 2); % number of neurons in target region

% append index to firing rate table
source_FR = [source_FR; 1:numsource]; % append index to firing rate table
target_FR = [target_FR; 1:numtarget];
FR_column = 1; index_column = 2;

% sort the firing rate (we've appended the index to the firing rate table
% which helps us trace back to the neurons)
firing_rate_ranking_source = sortrows(source_FR');
firing_rate_ranking_target = sortrows(target_FR');

% we transposed before sorting, so we transpose back
target_FR_ranked = firing_rate_ranking_target';
source_FR_ranked = firing_rate_ranking_source';

%% (2) create hist-counts of bins that capture FR distributions discretely
% initialize neuron counts, pattern matrices, and index trackers
nPattern = numel(X_target);
Y_source  = cell(size(X_target));
Y_target  = cell([2, size(X_target)]);

% officially binning the firing rates into nBins
histcounts_target = histcounts(firing_rate_ranking_target(:,1), nBins);
histcounts_source = histcounts(firing_rate_ranking_source(:,1), nBins);

% ESTABLISH MAX CELLS PER BIN
% initialize the bin counts and the sourceSmaller_marker
cellsToPick = zeros(1,nBins);
areaThatBounded_marker = zeros(1,nBins);
% iterate over the bins to pick indices and neurons from the two regions
for i = 1:nBins
    % Look at both bins
    target_count = histcounts_target(i);
    source_count = histcounts_source(i);
    % Pick the MAX CELLS we will sample from the brain area with the least
    nCells = min(target_count, source_count);
    cellsToPick(i) = nCells;
    % Which area are we going to randomly sample?
    if target_count > source_count
        areaThatBounded_marker(i) = 1; % same source for this bin was smaller
    end
end


%% (3) iterate over the bins to pick indices and neurons from the two regions

% indexTargetCount = cumsum(cellsToPick);
nTarget = 0;
index_target = zeros(2,nTarget);

sourceCountTracker = 0;
targetCountTracker = 0;

sourceTracker = 0;
targetTracker = 0;


for i = 1:nBins
    
    target_count = histcounts_target(i);
    source_count = histcounts_source(i);
    % pick the smaller neuron per bin size
    pickCells_i = cellsToPick(i);
    if Opt.verbose
        disp("----------------------------------")
        disp("i " + i)
        disp("pickCells_i " + pickCells_i)
        disp("target tracker " + targetCountTracker)
        disp("----------------------------------")
    end
        
    % Source has fewer neurons in this bin
    if areaThatBounded_marker(i) == 0
        % This condition checks if the source region has fewer neurons in the
        % current bin. In this case, all neurons from this bin in the target
        % region are taken, and an equal number of neurons are randomly
        % selected from the source region.
        indices_target = ... 
            target_FR_ranked(index_column,... 
                           targetTracker+1:targetTracker+pickCells_i);
        
        % randomly select the same # of cells from the source population
        random_ordering = randperm(source_count);
        indices_source = source_FR_ranked(index_column,... 
            random_ordering(1:pickCells_i) + sourceTracker);
        
        for j = 1:nPattern
            Y_target{1,j} = [Y_target{1,j}; X_source{j}(indices_source,:, :)];
            Y_target{2,j} = [Y_target{2,j}; X_target{j}(indices_target,:, :)];
        end
        
    % Target has fewer neurons in this bin
    else
        % This is the case when the target region has fewer neurons in the
        % current bin. In this case, all neurons from this bin in the source
        % region are taken, and an equal number of neurons are randomly
        % selected from the target region.
        indices_source = source_FR_ranked(index_column,...
            sourceTracker+1:sourceTracker+pickCells_i);
        
        % randomly select the same # of cells from the source population
        %         disp(indices_source + " same source less")
        % index the indices of source population first
        random_ordering = randperm(target_count);
        indices_target = target_FR_ranked(index_column,... 
            random_ordering(1:pickCells_i) + targetTracker);
        % select neurons from the source region
        for j = 1:nPattern
            Y_target{1,j} = [Y_target{1,j}; X_source{j}(indices_source,:, :)];
            Y_target{2,j} = [Y_target{2,j}; X_target{j}(indices_target,:, :)];
        end
        
    end
    
    assert(targetCountTracker <= numtarget, "targetCountTracker out of bounds")
    if  any(indices_target > numtarget)
        disp("An index is out of bounds")
        keyboard
    end

    % update the index trackers
    index_target(2,targetCountTracker+1: targetCountTracker+pickCells_i) = indices_target;
    targetCountTracker = targetCountTracker + pickCells_i;
    
    % update the index trackers
    index_target(1,sourceCountTracker+1: sourceCountTracker+pickCells_i) = indices_source;
    sourceCountTracker = sourceCountTracker + pickCells_i;
    
    targetTracker = targetTracker + target_count;
    sourceTracker = sourceTracker + source_count;
    
    nTarget = nTarget + pickCells_i;
    if numsource - nTarget <= numsource/3
        break
    end

end

%% finally, taking care of the source
nSource = numsource-nTarget;
index_source = setdiff((1:numsource),index_target(1,:));
for i = 1:nPattern
    Y_source{i} = X_source{i}(index_source,:, :);
end
