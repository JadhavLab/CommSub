function X = fieldSet(X, field, values, varargin)
% Set field values in a struct array
% Opposite of fieldGet()
%
% Usage:
%   X = nd.fieldSet(X, field, values, varargin)
%
% Inputs:
%   X: struct array to set fields in
%   field: field name to set
%   values: values to set field to
%       if values is a struct array, then assign X.(field) from values.(field)
%       if values is a cell array, then works like a set of calls to fieldSet()
%   varargin: passed to fieldGet()

assert(isstruct(X), 'X must be a struct');
if iscell(field) || ...
        (isstring(values) && numel(values) > 1)
    % if field is a cell array, then values must be a cell array of the same
    % size, and we call fieldSet() for each pair
    for i = 1:numel(field)
        X = nd.fieldSet(X, field{i}, values{i}, varargin{:});
    end
    return;
elseif isstruct(values)
    % if values is a struct array, then assign X.(field) from values.(field)
    % if values is a struct array, but doesn't have the field, then assign
    % the entire struct at i to the field
    if isfield(values, field)
        for i = 1:numel(X)
            X(i).(field) = values(i).(field);
        end
    else
        % else if values is a struct array, but doesn't have the field, then
        % assign the entire struct at i to the field
        for i = 1:numel(X)
            X(i).(field) = values(i);
        end
    end
else
    % else values is a scalar, so assign it to all elements of X
    for i = 1:numel(X)
        X(i).(field) = values(i);
    end
end
