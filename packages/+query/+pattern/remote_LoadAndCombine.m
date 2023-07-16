function [Patterns, otherData] = tmpLoadAndCombine(keys, varargin)
% [Patterns, otherData] = tmpLoadAndCombine(keys)
%
% Server analogue off loadAndCombine() - load and combine patterns from
% server into a single Nd array of patterns.
%
% INPUTS
%   keys    1xN cell array of keys to load
%
% OUTPUTS
%   Patterns    Nd array of patterns

ip = inputParser();
ip.addParameter('verbose', [], @(x) ischar(x) || isstring(x));
ip.parse(varargin{:});
Opts = ip.Results;

% rsycn files from server to local
kCount = 0;
mCount = 0; % number of missing files
for key = progress(keys(:)','Title','Remote->Local')
    kCount = kCount + 1;
    serverFile = fullfile(serverdefine, 'hash', key + ".mat");
    serverFile = serverFile.replace('//','/');
    localFile  = fullfile('/','tmp', key+".mat");
    disp("Downloading " + serverFile + " to " + localFile);
    command = "rsync --partial --progress -avu " + serverFile + " " + localFile;
    disp(command);
    system(command);
    if ~exist(localFile,'file')
        mCount = mCount + 1;
    end
end
disp("Missing " + mCount + " files");

% load and combine
otherData = cell(1,numel(keys));
kCount = 0;
P = cell(1,numel(keys));
for key = progress(keys(:)','Title','Loading keys')
    kCount = kCount + 1;
        keyboard
    tmp = load(localFile);
    P{kCount} = tmp.Patterns;
    tmp = rmfield(tmp,'Patterns');
    otherData{kCount} = tmp;
    if ~isempty(Opt.verbose) &&  isfield(tmp,'Option')
        t=struct2table(tmp.Option);
        disp("Option for " + key)
        disp("----------------")
        disp(t(:, Opt.verbose))
        disp("----------------")
    end
end

Patterns = ndb.toNd(squeeze(P));
