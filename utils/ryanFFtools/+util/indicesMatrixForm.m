function D = indicesMatrixForm(sz)

d = arrayfun(@(x) 1:x, sz, 'UniformOutput' , false);
D = cell(size(d));
[D{:}] = ndgrid(d{:});
D = cellfun(@(x) x(:), D, 'UniformOutput', false);
D = cat(2, D{:});
