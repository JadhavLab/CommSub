function [source_same, source_diff, target, same_index, different_source, t_index]...
                      = threeGroupsSplit(stFR, dtFR, X_same, X_different, nBins)
%
%
% this function splits the neuronal population into three groups, two
% sources and one target
%
% st stands for same as target
%
% dt stands for different from target
%
% it will divide the region serving as target into two and fr-match same
% number of source and target, and try its best to pick the same number of
% fr-matched neurons from the other region different than the target
%
% sort the firing pattern matrices by firing rate so as to index neurons


numSt = size(stFR,2);
numDt = size(dtFR,2);

stFR = [stFR; 1:numSt];
dtFR = [dtFR; 1:numDt];

fr_ranking_same = sortrows(stFR');
fr_ranking_different = sortrows(dtFR');

histcounts_same = histcounts(fr_ranking_same(:,1), nBins);
histcounts_different = histcounts(fr_ranking_different(:,1), nBins);

%% Iterate over the bins to pick indices and neurons from the two regions
stTracker = 0;
dtTracker = 0;
same_index = [];
different_source = [];
t_index = [];

% Iterate each of the firing rate bins
for i = 1:nBins
    same_count = histcounts_same(i);
    dt_count = histcounts_different(i);
    
    same_index_all = fr_ranking_same(stTracker+1: stTracker+same_count, 2);
    different_source_all = fr_ranking_different(dtTracker+1: dtTracker+dt_count, 2);
    rand_same = same_index_all(randperm(same_count));
    rand_different = different_source_all(randperm(dt_count));
    
    same_index =  [same_index, rand_same(1:ceil(same_count/2))'];
    t_index =  [t_index, rand_same(ceil(same_count/2)+1 : same_count)'];
    
    if dt_count>= same_count/2
        different_source = [different_source, rand_different(1:ceil(same_count/2))'];
    else
        different_source = [different_source, rand_different(1:ceil(dt_count))'];
    end
    
    stTracker = stTracker+same_count;
    dtTracker = dtTracker+dt_count;
end

%% Map back to the pattern matrices
nPatterns = size(X_same,2);

source_same = cell(1,nPatterns);
source_diff = cell(1,nPatterns);
target = cell(1,nPatterns);

for i = 1:nPatterns
    source_same{i} = X_same{i}(same_index,:);
    source_diff{i} = X_different{i}(different_source,:);
    target{i} = X_same{i}(t_index,:);
end

