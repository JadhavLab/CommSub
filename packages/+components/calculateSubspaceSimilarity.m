function [subspace, subspace_belonging] = calculateSubspaceSimilarity...
    (optDim, B_, spikeRaster, celllookup, sourceRegion, targetRegion, varargin)

ip = inputParser;
ip.addParameter('source_index', []);
ip.addParameter('target_index', []);

ip.parse(varargin{:});
opt = ip.Results;

% outputs which subspace overlaps the most across time for the spike
% ratser matrix

% stitch together the spiking matrices by the ranks
[u,~,v] = svd(B_);
subspace = [];


for iRank = 1:optDim
    curr_source = u(:,iRank);
    curr_target = v(:,iRank);
    
    curr_cells = [curr_source', curr_target'];
    subspace = [subspace; curr_cells];
end



%%
if ~isempty(opt.source_index) | ~isempty(opt.target_index)
    % get the speicifc cells from the whole spike raster by the source/target
    % indices

    selected_spikeRaster = [];
    for i = 1:numel(opt.source_index)
        for j = 1:height(celllookup)
            if table2array(celllookup(j,1)) == sourceRegion && table2array(celllookup(j,3)) == opt.source_index(i)
                selected_spikeRaster = [selected_spikeRaster; spikeRaster(j,:)];
            end
        end
    end
    %%
    % why not selecting all the target?
    % also need to generate celllookup and store it
    for i = 1:numel(opt.target_index)
        for j = 1:height(celllookup)
            if table2array(celllookup(j,1)) == targetRegion && table2array(celllookup(j,3)) == opt.target_index(i)
                disp(opt.target_index(i))
                selected_spikeRaster = [selected_spikeRaster; spikeRaster(j,:)];
            end
        end
    end
    %% take the subspaces
    try
    subspace_belonging = subspace * selected_spikeRaster;
    catch
        keyboard
    end
else
    subspace_belonging = subspace * spikeRaster;
end

end

