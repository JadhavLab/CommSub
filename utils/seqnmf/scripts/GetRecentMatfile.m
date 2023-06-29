% Enter that folder and load up the data
% --------------------------------------
% Input
% ---------
% fieldstr, kws, timescale, rootFolder
% Output
% ------
% D

rootedfolder = @(x) fullfile(rootFolder, x);
folder = seqnmf_folder(fieldstr, kws{:}, 'timescale', timescale, 'usedate', false, 'seqStyle', '', 'epoch_type', epoch_type);
D = dir(rootFolder);
candidates = arrayfun(@(x) contains(D(x).name, folder), 1:numel(D));
D = D(candidates);
if exist('ParamCandidates','var')
    ParamCandidates{end+1} = D;
end
notFoundFlag = isempty(D);
D
if ~isempty(D)
    [~,q] =sort(string({D.name}), 'descend');
    D = D(q(1));
    folderfull = fullfile(D.folder, D.name);
    disp(['File found = ' char(fullfile(folderfull,seqStyle,'master_seqnmf'))]);
    Mfile = fullfile(folderfull,seqStyle,'master_seqnmf');
    clear D;
else
    Mfile = "";
end
if exist(Mfile + '.mat', 'file')
    M = matfile(Mfile); 
    if ~isfield(M,'data') 
        warning('... file is missing data')
    end
    if ~ismember('behavior', fieldnames(M))
        warning('... file is missing behavior')
    end
    missingFlag = ~ismember('behavior', fieldnames(M)) || ~ismember('data',fieldnames(M));
else
    folderfull = [];
    disp(['File NOT found = ' fullfile(seqStyle,'master_seqnmf')]);
    M = [];
    missingFlag = false;
end
%try; clear Mfile; end;
clear rootedfolder candidates
