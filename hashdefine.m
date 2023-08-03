function path = hashdefine(varargin)

if isunix
    path = '/Volumes/Ark/commsubspace/hash/';
end

if ~isempty(varargin)
    path = fullfile(path,varargin{:});
end
