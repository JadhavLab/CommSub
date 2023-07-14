function [varargout] = combinePatterns(keys, varargin)
% Takes a list of keys and loads into combined pattern across keys
%
% final shape : keyCount x partitions x direction x patterns

disp('Adding multiEpoch to path');
addpath('hash_FRfixed');


ip = inputParser;
ip.addParameter('fields', ["Patterns", "Behaviors", "Option"]);
ip.addParameter('removeEmpty', false);
ip.addParameter('as', 'struct'); % cell or struct
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
if strcmp(Opt.as,'struct')
    for field = fields'
        fcnt = fcnt + 1;
        varargout{1}.(field) = {};
    end
elseif strcmp(Opt.as,'cell')
    for field = fields'
        fcnt = fcnt + 1;
        varargout{fcnt} = {};
    end
end

for key = progress(keys','Title','Loading keys')

    kCount = kCount + 1;
    if ~endsWith(key,".mat")
        key = key + ".mat";
    end

    if ~isempty(Opt.RunsSummary)
        tbl = Opt.RunsSummary(Opt.RunsSummary.hash == replace(key,".mat",""),:);
        disp("Loading " + key + " from " + tbl.animal + ...
            " with hash " + tbl.hash + " and " + tbl.generateH + " patterns" ...
            + " <- " + tbl.timestamp);
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
            disp("Loading " + field + " from " + key);
            if strcmp(Opt.as,'struct')
                varargout{1}.(field){kCount} = tmp.(field);
            elseif strcmp(Opt.as,'cell')
                varargout{fcnt}{kCount} = tmp.(field);
            end
        end
    end

end
%keyboard

% 
fcnt = 0;
for field = fields'
    fcnt = fcnt + 1;
    if strcmp(Opt.as,'cell')
        if Opt.removeEmpty
            empty = cellfun(@isempty, varargout{fcnt});
            varargout{fcnt}(empty) = [];
        end
        if ~isempty(varargout{fcnt})
             varargout{fcnt} = ndb.toNd(squeeze(varargout{fcnt}));
        end
    elseif strcmp(Opt.as, 'struct')
        if Opt.removeEmpty
            empty = cellfun(@isempty, varargout{1}.(field));
            varargout{1}.(field)(empty) = [];
        end
        if ~isempty(varargout{1}.(field))
            varargout{1}.(field) = ndb.toNd(squeeze(varargout{1}.(field)));
        end
    end
end






