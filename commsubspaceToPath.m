if ~exist('codedefine_to_path__', 'var') || ~codedefine_to_path__
    addpath(genpath(fullfile(codedefine())));
    cd(codedefine());
end
codedefine_to_path__ = true;
