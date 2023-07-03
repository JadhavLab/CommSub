function Dat = reshapeByLabels(Dat, dim, labels)
% Supposed we have an a x b x c x d structure, and
% on of those dimensions dim is labeled by labels.
% This function reshapes the data to be
% ... x (dim / nUniqueLabels) x ...

% Get the number of unique labels
[uLabels, ia, ic] = unique(labels);

% Assert that the labels have equal counts in ic
assert(all(histcounts(ic, 1:max(ic)+1) == length(ic) / length(uLabels)));
% Assert that all the labels increment by 1
dia = diff(ia);
assert(all(unique(dia) == dia), ...
    'Labels must equally spaced and equal in number');

% Reshape the data
Dat = reshape(Dat, [size(Dat, 1:dim-1), ...
    length(uLabels), size(Dat,dim)/length(uLabels), ... expanded dims
    size(Dat, dim+1:ndims(Dat))]);

