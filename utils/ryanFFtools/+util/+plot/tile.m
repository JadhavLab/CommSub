function tile(func, tiledim, varargin)
% Executes a plot function over tiles
% function tile(func, tiledim, varargin)

stringLoc = cellfun(@(x) isstring(x) || ischar(x), varargin);
if any(stringLoc)
    V = varargin(stringLoc:end);
    varargin(stringLoc:end) = [];
else
    V = {};
end
ip = inputParser;
ip.addParameter('kws', {});
ip.parse(V{:})
Opt = ip.Results;

for d = progress(util.indicesMatrixForm(size(varargin{1}, tiledim))')

    ind = repmat(':', 1, ndims(varargin{1}));
    ind = num2cell(ind);
    ind{tiledim} = d;

    d = num2cell(d);
    nexttile(d{:});
    V = cellfun(@(x) squeeze(x(ind{:})), varargin, 'UniformOutput',false);
    if ~iscell(func) && ~isstruct(func)
        func(V{:}, Opt.kws{:});
    elseif iscell(func)
    elseif isstruct(func)
    end

end
