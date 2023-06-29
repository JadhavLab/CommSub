function [inds, logical_inds] = field(data, field, varargin)
% INC SIEVE.FIELD produces the indices of the axes that cut the proper field out of the data.data structure (or the data.W* structures.

% Use the labeler function to find the indices that should correspond to our field of interest
ip = inputParser;
ip.addParameter('removefield',{});
ip.parse(varargin{:})
opt = ip.Results;

[~,~,pts,fieldPts] = labeler.wholeLabel(data, data.fields(~ismember(data.fields,opt.removefield))); % Get the name of all the points of the property axis of data.data without downsampling
logical_inds = fieldPts == string(field);
inds = find(logical_inds);
