function x = figuredefine(varargin)
% x = figuredefine(varargin)
%
%   Defines the path to a figure folder.  The path is relative to the codedefine
%   folder.  The varargin inputs are the subfolders and the filename.  The
%   subfolders and filename are optional.
%
% other options:
%   -permfolder:  sets the permanent subfolder
%   -clearpermfolder:  clears the permanent subfolder
%   -getpermfolder:  returns the permanent subfolder
%   -creation:  if true, figuredefine always creates the folder if it doesn't exist
%

persistent permfolder
if varargin{1} == "-permfolder"
    permfolder = varargin(2:end);
    return
elseif varargin{1} == "-clearpermfolder"
    permfolder = [];
    return
elseif varargin{1} == "-getpermfolder"
    x = permfolder;
    return
end
persistent creation
if ~exist('creation','var')
    creation = false;
end
if varargin{1} == "-creation"
    creation = varargin{2};
    return
end

if ~exist('permfolder','var') || isempty(permfolder)
    x = string(fullfile(codedefine,'figures', varargin{:}));
else
    x = string(fullfile(codedefine,'figures', permfolder{:}, varargin{:}));
end

if creation
    [~,base,ext] = fileparts(x);
    if ~exist(x,'dir') && ~contains(base,'.')
        mkdir(x)
    end
end
