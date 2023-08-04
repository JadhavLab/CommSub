function [subspaceDist, rowVar] = load(k, normalize, measurement, selectGenH)

figureFolder = figuredefine("subspaceAngle","type="+selectGenH);
analysisName = sprintf('subspaceDist_K=%d_norm=%d_measure=%s.mat',...
    k, normalize, measurement);
if ~exist(figuredefine("subspaceAngle","type="+selectGenH), 'dir')
    mkdir(figuredefine("subspaceAngle","type="+selectGenH))
end
analysisName = figuredefine("subspaceAngle","type="+selectGenH,...
                            analysisName);
load(analysisName, 'subspaceDist', 'rowVar');
% ----------------------------
% Normalize and plot distances
% ----------------------------
% f = figc("subpsace distances, " + analysisName);
for i = 1:size(subspaceDist,1)
    for j = 1:size(subspaceDist,2)
        % keyboard
        if i==j
            mu1 = subspaceDist(i,:);
            mu2 = subspaceDist(:,j);
            subspaceDist(i,j) = mean([mu1(:) mu2(:)], 'all','omitnan');
        end
    end
end
