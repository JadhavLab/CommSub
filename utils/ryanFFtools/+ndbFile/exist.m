function E = exist(animID, datatype, level)
% True if a file with the animID, datatype, level exists

if nargin < 3
    level = [];
end

folder = ndbFile.folder(animID, datatype);

if isempty(level)
    E = ~isempty(dir(fullfile(folder, sprintf('%s%s*',animID,datatype))));
else
    error("Level not implemented yet");
end

