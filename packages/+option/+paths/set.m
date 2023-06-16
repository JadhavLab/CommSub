function set()
% option.paths.set() sets paths for a given machine

% TODO Move this paragraph into codedefine()
if ispc
    %pass
elseif ismac
    %pass
else
    addpath('/Volumes/MATLAB-Drive/Shared/');
end

persistent pathSet
if ~exist('pathSet','var') || isempty(pathSet)
    pathSet = false;
end

if pathSet == true
    disp("Path already set!")
else
    disp("Setting code path")
    paths = codedefine();
    controlled_path_set(paths);
    disp("Setting data path")
    paths = datadefine();
    controlled_path_set(paths);
    pathSet = true;
end

function controlled_path_set(paths)

    for path = string(paths)
        for gp = genpath(x)
            try
                rmpath(gp) % in case already in path
            end
            addpath(gp)
        end
     end
