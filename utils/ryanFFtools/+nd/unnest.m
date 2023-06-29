function [obj, fields] = unnest(obj, nestfield, fields)
% Unnests all fields in nestfield, or optionally, certain fields

if nargin < 3
    fields = [];
end

if numel(obj) > 1
    obj = arrayfun(@(x) nd.unnest(x, nestfield, fields), obj, 'UniformOutput', false);
    obj = ndb.toNd(obj);
else
    if ~isfield(obj, nestfield) || ~isstruct(obj.(nestfield))
        obj = [];
        fields = [];
        return;
    end

    nFields = numel(fieldnames(obj.(nestfield)));
    if nargin <= 2 || isempty(fields)
        fields = string(fieldnames(obj.(nestfield)));
    else
        fields = string(fields);
    end

    if numel(fields) == nFields
        deleteNestField  = true;
    else
        deleteNestField = false;
    end

    for field = fields(:)'
        if isfield(obj.(nestfield), field)
            obj.(field) = obj.(nestfield).(field);
            obj.(nestfield) = rmfield(obj.(nestfield),field);
        end
    end

    if deleteNestField
        obj = rmfield(obj, nestfield);
    end
end
