function [s_hpc, s_pfc, t_hpc, t_pfc,s_hpc_index,s_pfc_index...
    ,t_hpc_index,t_pfc_index] = twoWaySplit(pfcFR, hpcFR, X_hpc, X_pfc, nBins)
% this function splits the neuronal population into two ways, by dividing the
% number of neurons equally matched by firing rate.

numhpc = size(hpcFR,2);
numpfc = size(pfcFR,2);
s_hpc_index = [];
s_pfc_index = [];
t_hpc_index = [];
t_pfc_index = []; 

pfcFR = [pfcFR; 1:numpfc];
hpcFR = [hpcFR; 1:numhpc];

fr_ranking_pfc = sortrows(pfcFR');
fr_ranking_hpc = sortrows(hpcFR');

% get the firing rate of pfc and hpc neurons divided into bin numbers
pfc_firing = histcounts(fr_ranking_pfc(:,1), nBins);
hpc_firing = histcounts(fr_ranking_hpc(:,1), nBins);

%% select neurons from both regions ranked by the firing rate bins
hpcTracker = 0;
pfcTracker = 0;

for i = 1:nBins
    hpc_count = hpc_firing(i);
    pfc_count = pfc_firing(i);
    
    % randomize the indices
    hpc_index = fr_ranking_hpc(hpcTracker+1:hpcTracker+hpc_count,2);
    rand_hpc = hpc_index(randperm(hpc_count));
    pfc_index = fr_ranking_pfc(pfcTracker+1:pfcTracker+pfc_count,2);
    rand_pfc = pfc_index(randperm(pfc_count));
    
    % append to respective groups
    s_hpc_index =  [s_hpc_index, rand_hpc(1:ceil(hpc_count/2))'];
    t_hpc_index =  [t_hpc_index, rand_hpc(ceil(hpc_count/2)+1 : hpc_count)'];
    s_pfc_index =  [s_pfc_index, rand_pfc(1:ceil(pfc_count/2))'];
    t_pfc_index =  [t_pfc_index, rand_pfc(ceil(pfc_count/2)+1 : pfc_count)'];
    
    % increment tracker
    hpcTracker = hpcTracker+hpc_count;
    pfcTracker = pfcTracker+pfc_count;

end

%% map back to the pattern matrices

numResult = size(X_hpc,2);
s_hpc = cell(1,numResult);
s_pfc = cell(1,numResult);
t_hpc = cell(1,numResult);
t_pfc = cell(1,numResult);

for i = 1:numResult
    s_hpc{i} = X_hpc{i}(s_hpc_index,:);
    s_pfc{i} = X_pfc{i}(s_pfc_index,:);
    
    t_hpc{i} = X_hpc{i}(t_hpc_index,:);
    t_pfc{i} = X_pfc{i}(t_pfc_index,:);
end

end

