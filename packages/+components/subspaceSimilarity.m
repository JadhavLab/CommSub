function [total_subspaces, cell_subspaces] = ...
        subspaceSimilarity_v2(dim, B_cell, ...
    spikeRateMatrix, celllookup, source, target, varargin)
% function [total_subspaces, cell_subspaces] = ...
%     subspaceSimilarity_v2(dim, B_cell, ...
%     spikeRateMatrix, celllookup, source, target, varargin)
% Calculates the subspace similarity for the rate matrix to a given set of 
% subpsace components
% INPUTS:
%   dim: dimensionality of the subspace
%   B_cell: cell array of the subspace components
%   spikeRateMatrix: matrix of spike rates
%   celllookup: cell lookup table
%   source: source region
%   target: target region
%   varargin: optional arguments
%       source_index: index of the source region
%       target_index: index of the target region
% OUTPUTS:
%   total_subspaces: combined subspace
%   cell_subspaces: cell array of the subspace for each cell


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
if  ~isempty(opt.source_index) || ~isempty(opt.target_index)
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
