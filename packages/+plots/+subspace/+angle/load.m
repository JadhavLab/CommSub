function [subspaceDist, frobDist, rowVar] = load(k, normalize, measurement, selectGenH)

figureFolder = figuredefine("subspaceAngle","type="+selectGenH);
analysisName = sprintf('subspaceDist_K=%d_norm=%d_measure=%s.mat',...
    k, normalize, measurement);
if ~exist(figuredefine("subspaceAngle","type="+selectGenH), 'dir')
    mkdir(figuredefine("subspaceAngle","type="+selectGenH))
end
analysisName = figuredefine("subspaceAngle","type="+selectGenH,...
                            analysisName);
disp("loading " + analysisName)
load(analysisName, 'subspaceDist', 'frobDist', 'rowVar');
% ----------------------------
% Normalize and plot distances
% ----------------------------
% f = figc("subpsace distances, " + analysisName);
for k = 1:size(subspaceDist,1)
    for i = size(subspaceDist,2)
    for j = size(subspaceDist,3)
        subspaceDist(k,i,j) = mean([subspaceDist(k, i, j), subspaceDist(k, j, i)]);
    end
    end
end

