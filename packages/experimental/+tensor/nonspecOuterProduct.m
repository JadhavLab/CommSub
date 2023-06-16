function Z = nonspecOuterProduct(X, varargin)
%
% creates a tensor outer product of two area-specific tensors, where time slices
% of a given output tensor are the outer product pieces that would be averaged
% to obtain a B regression matrix.

Z = tensor.areawiseOuterProduct(X, X, varargin{:});
