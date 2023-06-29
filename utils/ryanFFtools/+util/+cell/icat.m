function Tnew = icat(Tcell, varargin)
% T = function icat(Tcell)
% A field insensitive version of cat() that's immune to differences in
% field number --- for cating cells of table or cells of struct
%
% Method will either concatonate only the fields that match for all
% tables (intersect) or concatonate all fields, including unmatched
% (union), where tables that lack a field have those missing fields as
% nan.

dim = [];
if ~isempty(varargin) && isnumeric(varargin{1})
    dim = varargin{1};
    varargin(1) = [];
end

ip = inputParser;
ip.addParameter('fieldCombine', 'intersect');
ip.addParameter('removeEmpty', false);
ip.parse(varargin{:})
Opt = ip.Results;


if isempty(dim)
    if Opt.removeEmpty
        Tcell = util.cell.removeEmpty(Tcell);
    end
    fields = util.cell.fieldnames(Tcell, Opt.fieldCombine);
    Tcell = cellfun(@(x) util.table.select(x, fields), Tcell, 'UniformOutput', false);
    Tnew = cat(1, Tcell{:});
else
    inds = nd.indicesMatrixForm(Tcell);
    inds = unique(inds(:, setdiff(1:ndims(inds), dim)), 'rows');
    inds = num2cell(inds);
    insert = @(x, n) cat(2, x(:,1:(n-1)), repmat({':'}, size(x,1),1), x(:,n:end));
    inds = insert(inds, dim);
    indsNew = inds;
    indsNew(:,dim) = [];

    Tnew = {};
    for ii = progress(1:size(inds,1), 'Title', sprintf('Concat over dim=%d', dim))
        ind     = inds(ii,:);
        indNew = indsNew(ii,:);
        Tnew{indNew{:}} = util.cell.icat(Tcell(ind{:}));
        [Tcell{ind{:}}] = deal([]);
    end
end
