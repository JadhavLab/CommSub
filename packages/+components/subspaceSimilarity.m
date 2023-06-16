function [total_subspaces, cell_subspaces] = subspaceSimilarity_v2(dim, B_cell, ...
    spikeRateMatrix, celllookup, source, target, varargin)
% Gives the cell subspace and combined subspace
ip = inputParser;
ip.addParameter('source_index', []);
ip.addParameter('target_index', []);

ip.parse(varargin{:});
opt = ip.Results;



if nargin < 5
    region = "CA1";
end

subspaces = {};
if  ~isempty(opt.source_index) | ~isempty(opt.target_index)
    for i = 1:numel(B_cell)
        [~, cell_subspaces{i}] = components.calculateSubspaceSimilarity...
            (dim, B_cell{i}, spikeRateMatrix, celllookup,source,target,...
            'source_index', opt.source_index,...
            'target_index', opt.target_index);
    end
else
    for i = 1:numel(B_cell)
        [~, cell_subspaces{i}] = components.calculateSubspaceSimilarity...
            (dim, B_cell{i}, spikeRateMatrix, celllookup,"CA1");
    end
end
total_subspaces = cat(1, cell_subspaces{:});
