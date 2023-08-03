function save(subspaceDist, rowVar, k, normalize, measurement, selectGenH)

analysisName = sprintf('subspaceDist_K=%d_norm=%d_measure=%s.mat',...
    k, normalize, measurement);
if ~exist(figuredefine("subspaceAngle","type="+selectGenH), 'dir')
    mkdir(figuredefine("subspaceAngle","type="+selectGenH))
end
analysisName = figuredefine("subspaceAngle","type="+selectGenH, analysisName);
save(analysisName, 'subspaceDist', 'rowVar');
