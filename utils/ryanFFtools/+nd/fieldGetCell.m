function out = fieldGet(X, field)
% Gets a single field from the nd struct, and stores it into a cell in the same nd shape

out = cell(size(X));
indices = nd.indicesMatrixForm(X);
for index = indices'
    I = num2cell(index);
    out{ I{:} } = X( I{:} ).(field);
end
