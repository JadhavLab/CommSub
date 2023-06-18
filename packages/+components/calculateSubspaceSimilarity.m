function [subspace, subspace_belonging] = calculateSubspaceSimilarity...
    (optDim, B_, spikeRaster, celllookup, sourceRegion, targetRegion, varargin)
%CALCULATESUBSPACESIMILARITY calculate the similarity of the firing rate
% activity of the cells to the subspace provided in B_
%
% Inputs:
%   optDim - the number of dimensions to use for the subspace, Optimal Dimension
%   B_ - the subspace
%   spikeRaster - the spike raster matrix
%   celllookup - the cell lookup table
%   sourceRegion - the source region, string 
%   targetRegion - the target region, string
%   varargin - 'source_index' - the index of the source cells
%              'target_index' - the index of the target cells
%
% Outputs:
%   subspace - the subspace
%   subspace_belonging - the similarity of the firing rate activity of the
%   cells to the subspace

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


% build a row of the expected source and target neuron activity
for iRank = 1:optDim
    curr_source = u(:,iRank);
    curr_target = v(:,iRank);
    
    curr_cells = [curr_source', curr_target'];
    subspace = [subspace; curr_cells];
end



% calculate the similarity of the firing rate activity of the cells to the
% subspace
if ~isempty(opt.source_index) || ~isempty(opt.target_index)

    % grab the correct cells specified by the source and target index
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

