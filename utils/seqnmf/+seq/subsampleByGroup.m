function data = subsampleByGroup(data, fields, subsample, varargin)
% SUBSAMPLEBYGROUP subsamples some percentage within each group
%
% Inputs
% ------
% data : seqnmf struct
%   See calling code for information about how one would struture this struct
% fields : cell
%   List of fields in the struct to apply the subsampling operation onto
% subsample : scalar
%   A scalar number encoding the percentage of samples to take.
%
% Rewrite
% -------
% No longer makes any calls to any subfunctions.
% 
% Returns
% -------
% data : struct
%   Afformentioned data structure with subsampled data.

ip = inputParser;
ip.addParameter('storeOrig', false);
ip.parse(varargin{:});
opt = ip.Results;

data.subsample = sprintf('%2.1f percent', 100*subsample); 

if subsample == 1
    return
end

G = findgroups(data.animalcnt, data.epoch);
fields = [fields, {'t'}];
if ismember('trajdist', fields)
    fields = [fields, {'orig_trajdist'}];
end
if ismember('lindist', fields)
    fields = [fields, {'orig_lindist'}];
end
if ismember('aglindist', fields)
    fields = [fields, {'orig_aglindist'}];
end
if ismember('velocity', fields)
    fields = [fields, {'orig_velocity'}];
end
if ismember('reward', fields)
    fields = [fields, {'orig_reward'}];
end
if ismember('trajbound', fields)
    fields = [fields, {'orig_trajbound'}];
end
if isfield(data,'data')
    fields = ['data' fields];
end
fields

% --------- Subset the data ------------------
fields_emptyArray = arrayfun(@(x) [], 1:numel(fields), 'UniformOutput', false);
F = [fields; fields_emptyArray];
subData = struct(F{:});
I = zeros(size(data.t));
for group = unique(G)
    %group
    group_filter = G == group;
    t             = data.t(1, group_filter);
    subsampleSize = round( numel(t) * subsample );
    startpoint    = randi(numel(t) - subsampleSize, 1);
    endpoint      = startpoint + subsampleSize;
    for field = fields
        field = field{1};
        if isequal(field, 't')
            subData.t = [subData.t, ...
                         t(1, startpoint:endpoint)];
             I = I | (group_filter & ismember(data.t, t(1,startpoint:endpoint)));
        elseif isequal(field, 'data')
            d = data.data(:, group_filter);
            subData.data = [subData.data, ... 
                         d(:, startpoint:endpoint)];

        else
            d = data.(field)(group_filter, :);
            subData.(field) = [subData.(field); ...
                              d(startpoint:endpoint,:)];
            b = size(d, 1);
            assert(numel(t) == b, "Size mismatch")
        end
    end
end

data.subsample_logical = I;
data.subsample_indices = find(I);
%data.inds = data.subsample_indices;

for field = fields
    field = field{1};
    assert(any( size(data.(field)) > size(subData.(field)) ));
    data.(field) = single(subData.(field));
end
