function [sourceCellsInvolved, targetCellsInvolved] = ...
    getParticipatingCells(B_, index_source, index_target)

% B_ : the predicition matrix where entries of 0 correspond to cells that
% do not participate in a dimension of population interaction

% output the list of source and target cells involved in a certain
% dimension of population interaction in terms of their indices

%%
[nSource, nTarget] = size(B_);
zeroSource = [];
zeroTarget = [];

%% get all the participating source cells
for i = 1:nSource
    if sum(abs(B_(i,:)) > eps('double')) > 0
        zeroSource = [zeroSource, i];
    end
end

%% get all the participating target cells
for i = 1:nTarget
    if sum(abs(B_(:,i)) > eps('double')) > 0
        zeroTarget = [zeroTarget, i];
    end
end

%% map back to the index list
sourceCellsInvolved = index_source(zeroSource);
targetCellsInvolved = index_target(zeroTarget);
end

