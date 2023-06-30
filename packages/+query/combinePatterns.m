function [varargout] = combinePatterns(keys, varargin)
% Takes a list of keys and loads into combined pattern across keys
%
% final shape : keyCount x partitions x direction x patterns

disp('Adding multiEpoch to path');
addpath('hash_FRfixed');


ip = inputParser;
ip.addParameter('fields', ["Patterns", "Behaviors", "Option"]);
ip.addParameter('removeEmpty', false);
ip.addParameter('RunsSummary', table(), @istable);
ip.parse(varargin{:})
Opt = ip.Results;
fields = string(Opt.fields(:));

kCount = 0;
P = cell(1,numel(keys));
B = cell(1,numel(keys));

otherData = cell(1,numel(keys));
keys = keys(:);

fcnt = 0;
for field = fields'
    fcnt = fcnt + 1;
    varargout{fcnt} = {};
end

for key = progress(keys','Title','Loading keys')

    kCount = kCount + 1;
    if ~endsWith(key,".mat")
        key = key + ".mat";
    end

    if ~isempty(Opt.RunsSummary)
        tbl = Opt.RunsSummary(Opt.RunsSummary.hash == replace(key,".mat",""),:);
        disp("Loading " + key + " from " + tbl.animal + ...
            " with hash " + tbl.hash + " and " + tbl.generateH + " patterns");
    end
    if ~exist(key,'file')
        warning('File not found: %s', key);
        continue
    end

    tmp = matfile(key);

    fcnt = 0;
    for field  = fields'
        fcnt = fcnt + 1;
        if ismember(field, fieldnames(tmp))
            varargout{fcnt}{kCount} = tmp.(field);
        end
    end

end
%keyboard

% 
fcnt = 0;
for field = fields'
    fcnt = fcnt + 1;
    if Opt.removeEmpty
        empty = cellfun(@isempty, varargout{fcnt});
        varargout{fcnt}(empty) = [];
    end
    if ~isempty(varargout{fcnt})
         varargout{fcnt} = ndb.toNd(squeeze(varargout{fcnt}));
    end
end






