function varargout = loadAndCombine(keys, varargin)
% Takes a list of keys and loads into combined pattern across keys
%
% final shape : keyCount x partitions x direction x patterns

varargout = query.combinePatterns(keys, varargin{:});

% kCount = 0; % tracks the key number
% mCount = 0; % tracks missing keys
% P = cell(1,numel(keys));
% otherData = cell(1,numel(keys));
% for key = progress(keys(:)','Title','Loading keys')
%     kCount = kCount + 1;
%     file = fullfile(datadefine, 'hash', key + ".mat")
%     if exist(file,'file')
%         disp("Loading " + key)
%         if nargin > 1
%             disp(tbl(tbl.hash == key,["animal","date","generateH"]))
%         end
%         tmp = load(file);
%         P{kCount} = tmp.Patterns;
%         tmp = rmfield(tmp,'Patterns');
%         otherData{kCount} = tmp;
%     else
%         mCount = mCount + 1;
%         fprintf('Missing key %s\n', key);
%     end
% end
%
% Patterns = ndb.toNd(squeeze(P));
