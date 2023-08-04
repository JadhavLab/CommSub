function subspaceDist = getSubspaceDist(subspaceDist, measurement, normalize)
% GETSUBSPACEDIST calculates the distance between subspaces

if contains(measurement, "simil")
    subspaceDist = squeeze(median(subspaceDist,1));
    % subspaceDist = subspaceDist + abs(min(subspaceDist,[],'all'));
    nondiag = ~eye(size(subspaceDist,[2,3])); % added 2023, ry
    nondiag = repmat(nondiag,[size(subspaceDist,1),1,1]);
    maxVal = max(subspaceDist(nondiag),[],'all');
    subspaceDist = maxVal - subspaceDist;
elseif contains(measurement, "dissimil")
    subspaceDist = squeeze(median(subspaceDist,1));
else
    error("You spelled it wrong")
end
if normalize
    nondiag = ~eye(size(subspaceDist)); % added 2023, ry
    subspaceDist = (subspaceDist-min(subspaceDist(nondiag),[],'all')) ...
    ./(max(subspaceDist(nondiag),[],'all')-min(subspaceDist(nondiag),[],'all'));
else
    % NOTHING :)
end
