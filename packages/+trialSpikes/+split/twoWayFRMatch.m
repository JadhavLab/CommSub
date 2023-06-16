function [X_source, X_target, nSource, nTarget, index_source, index_target] =...
    matchFRinDiscreteRanges(X_samesourceRegion, X_diffsourceRegion, ...
    samesource_FR, diffsource_FR, nBins)
%
% this function matches the firing rates of two target region with discrete
% bins from the FR histogram for every pattern
% it will make the prediction matrix more populated by having two regions
% with neurons that have matching firing rate
%
%
% ------------------------------------------------------------------------%
% input: the firing pattern matrices from the two regions                 %
% output: the matched firing pattern matrices from the two target regions %
%         the number of source and target population eventually decided   %
%         and the neuronal indices of those used                          %
% ------------------------------------------------------------------------%
%

%% sort the firing pattern matrices by firing rate so as to index neurons
% samesource_FR = hpcFR;
% diffsource_FR = pfcFR;
numSameSource = size(samesource_FR,2);
numDiffSource = size(diffsource_FR,2);


% append index to firing rate table
samesource_FR = [samesource_FR; 1:numSameSource];
diffsource_FR = [diffsource_FR; 1:numDiffSource];

% sort the firing rate and keep the index so as to track when binning
firing_rate_ranking_samesource = sortrows(samesource_FR');
firing_rate_ranking_diffsource = sortrows(diffsource_FR');

diffsource_ranking = firing_rate_ranking_diffsource';
samesource_ranking = firing_rate_ranking_samesource';
%% create hist-counts of bins that capture FR distributions discretely

% initialize neuron counts, pattern matrices, and index trackers
nSource = 0;
nTarget = 0;
numResult = numel(X_diffsourceRegion);
X_source  = cell(size(X_diffsourceRegion));
X_target  = cell([2, size(X_diffsourceRegion)]);

histcounts_diffsource = histcounts(firing_rate_ranking_diffsource(:,1), nBins);
histcounts_samesource = histcounts(firing_rate_ranking_samesource(:,1), nBins);

binCounts = zeros(1,nBins);
sameSourceSmaller_marker = zeros(1,nBins);

for i = 1:nBins
    % Look at both bins
    diffsource_count = histcounts_diffsource(i);
    samesource_count = histcounts_samesource(i);
    % Pikc the cell count from the brain area with the least
    nCells = min(diffsource_count, samesource_count);
    binCounts(i) = nCells;
    
    % Which area are we going to randomly sample?
    if diffsource_count > samesource_count
        sameSourceSmaller_marker(i) = 1; % same source for this bin was smaller
    end
end

indexTargetCount = cumsum(binCounts);
nTarget = 0;
index_target = zeros(2,nTarget);

%% iterate over the bins to pick indices and neurons from the two regions

samesourceCountTracker = 0;
diffsourceCountTracker = 0;

samesourceRegionTracker = 0;
diffsourceRegionTracker = 0;


for i = 1:nBins
    
    diffsource_count = histcounts_diffsource(i);
    samesource_count = histcounts_samesource(i);
    % pick the smaller neuron per bin size
    nCells = binCounts(i);
    
    % NON-MATCHING AREA SELECTION (SOURCE ~= TARGET)
    if sameSourceSmaller_marker(i) == 0
        % take all cells from that bin from the population diff than
        % source
        

         indices_diffsource = diffsource_ranking(2,diffsourceRegionTracker+1:diffsourceRegionTracker+nCells);
        
        % randomly select the same # of cells from the source population
        possible_indices_of_indices = randperm(samesource_count);
        indices_samesource = samesource_ranking(2,possible_indices_of_indices(1:nCells) + samesourceRegionTracker);
        
        for j = 1:numResult
            X_target{2,j} = [X_target{2,j}; X_diffsourceRegion{j}(indices_diffsource,:, :)];
            X_target{1,j} = [X_target{1,j}; X_samesourceRegion{j}(indices_samesource,:, :)];
        end
        
        
    % MATCHING AREA SELECTION (SOURCE ~= TARGET)
    else
        % take all cells from that bin from the population same as
        % source
        
        indices_samesource = samesource_ranking(2,samesourceRegionTracker+1:samesourceRegionTracker+nCells);
        
        % randomly select the same # of cells from the source population
        %         disp(indices_samesource + " same source less")
        % index the indices of source population first
        possible_indices_of_indices = randperm(diffsource_count);
        indices_diffsource = diffsource_ranking(2,possible_indices_of_indices(1:nCells) + diffsourceRegionTracker);
        % select neurons from the source region
        for j = 1:numResult
            X_target{1,j} = [X_target{1,j}; X_samesourceRegion{j}(indices_samesource,:, :)];
            X_target{2,j} = [X_target{2,j}; X_diffsourceRegion{j}(indices_diffsource,:, :)];
        end
        
    end
    
    assert(diffsourceCountTracker <= numDiffSource)
    
    if  any(indices_diffsource > numDiffSource)
   
        keyboard
    end
    index_target(2,diffsourceCountTracker+1: diffsourceCountTracker+nCells) = indices_diffsource;
    diffsourceCountTracker = diffsourceCountTracker + nCells;
    
    index_target(1,samesourceCountTracker+1: samesourceCountTracker+nCells) = indices_samesource;
    samesourceCountTracker = samesourceCountTracker + nCells;
    
    diffsourceRegionTracker = diffsourceRegionTracker + diffsource_count;
    samesourceRegionTracker = samesourceRegionTracker + samesource_count;
    
    nTarget = nTarget + nCells;
    if numSameSource - nTarget <= numSameSource/3
        break
    end
end

%% finally, taking care of the source

nSource = numSameSource-nTarget;

index_source = setdiff((1:numSameSource),index_target(1,:));

for i = 1:numResult
    X_source{i} = X_samesourceRegion{i}(index_source,:, :);
end


