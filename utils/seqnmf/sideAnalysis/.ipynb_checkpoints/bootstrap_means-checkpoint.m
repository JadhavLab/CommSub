function [varargout] = bootstrap_means(varargin)

opt = struct('UseParallel',true);
for v = 1:numel(varargin)
    V=varargin{v};
    varargout{v} = bootstrp(1000, @mean, V, 'Options', opt);
end
